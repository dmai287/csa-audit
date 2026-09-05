"""Filtered range/null-space projectors and the null-space audit statistics.

P_lambda = (A^H A + lambda I)^{-1} A^H A,   Q_lambda = I - P_lambda.
Neither matrix is formed; each application solves one Hermitian
positive-definite system by conjugate gradients using only forward/adjoint
products. As lambda -> 0+, P_lambda -> A^+ A and Q_lambda -> I - A^+ A.
"""
from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
from scipy import ndimage

from csa.physics.operator import SenseOperator


def real_inner(a: np.ndarray, b: np.ndarray) -> float:
    """Real part of the complex inner product <a, b> = sum conj(a) b."""
    return float(np.vdot(a, b).real)


def cg_solve(apply: Callable[[np.ndarray], np.ndarray], b: np.ndarray,
             tol: float = 1e-10, maxiter: int = 500) -> Tuple[np.ndarray, dict]:
    """Conjugate gradients for a Hermitian positive-definite operator."""
    x = np.zeros_like(b)
    r = b.copy()
    p = r.copy()
    rs = real_inner(r, r)
    rs0 = np.sqrt(rs)
    info = {"iterations": 0, "converged": rs0 == 0.0, "relative_residual": 0.0 if rs0 == 0 else 1.0}
    if rs0 == 0.0:
        return x, info
    for k in range(1, maxiter + 1):
        Ap = apply(p)
        alpha = rs / real_inner(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rs_new = real_inner(r, r)
        info["iterations"] = k
        info["relative_residual"] = float(np.sqrt(rs_new) / rs0)
        if np.sqrt(rs_new) <= tol * rs0:
            info["converged"] = True
            break
        p = r + (rs_new / rs) * p
        rs = rs_new
    return x, info


def range_filter(op: SenseOperator, x: np.ndarray, lam: float, tol: float = 1e-10,
                 maxiter: int = 500) -> np.ndarray:
    """P_lambda x."""
    b = op.normal(x)
    z, _ = cg_solve(lambda v: op.normal(v) + lam * v, b, tol=tol, maxiter=maxiter)
    return z


def null_filter(op: SenseOperator, x: np.ndarray, lam: float, tol: float = 1e-10,
                maxiter: int = 500) -> np.ndarray:
    """Q_lambda x = x - P_lambda x."""
    return x - range_filter(op, x, lam, tol=tol, maxiter=maxiter)


def noise_lambda(noise_sigma: float, X: np.ndarray) -> float:
    """lambda = sigma_eta^2 / var(|X|): the paper's noise-set rule for the filter."""
    return float(noise_sigma ** 2 / np.var(np.abs(X)))


def measurement_gain(op: SenseOperator, ell: np.ndarray) -> float:
    """kappa(ell) = ||A ell|| / ||ell||."""
    return float(np.linalg.norm(op.forward(ell)) / np.linalg.norm(ell))


def measurement_fraction(op: SenseOperator, ell: np.ndarray, lam: float, **cg) -> float:
    """mu_lambda(ell) = ||P_lambda ell||^2 / ||ell||^2 (model-independent)."""
    p = range_filter(op, ell, lam, **cg)
    return float(np.vdot(p, p).real / np.vdot(ell, ell).real)


def transfer_statistics(op: SenseOperator, dX: np.ndarray, ell: np.ndarray, lam: float,
                        **cg) -> Tuple[float, float]:
    """(t_R, t_N): measured and unmeasured transfer of the lesion through the model.

    t_R = <A dX, A ell> / ||A ell||^2 compares the arms in the measurement
    domain: it is exactly 1 for any reconstruction that enforces data
    consistency, whether or not the lesion survives, and it is what a residual
    gate checks. t_N = <Q dX, Q ell> / ||Q ell||^2 measures the part the prior
    had to supply; it is the statistic the residual cannot see.
    """
    if not np.any(ell):
        raise ValueError("lesion perturbation has zero energy")
    Ad, Al = op.forward(dX), op.forward(ell)
    t_R = real_inner(Ad, Al) / real_inner(Al, Al)
    Qd = null_filter(op, dX, lam, **cg)
    Ql = null_filter(op, ell, lam, **cg)
    t_N = real_inner(Qd, Ql) / real_inner(Ql, Ql)
    return float(t_R), float(t_N)


def false_structure_count(null_error: np.ndarray, background: np.ndarray, z_thr: float,
                          min_area_px: int, region: np.ndarray = None) -> int:
    """Connected components of |Q(err)| above mean + z_thr*sd of background, of size >= min_area."""
    a = np.abs(null_error)
    mu, sd = a[background].mean(), a[background].std()
    cand = a > mu + z_thr * sd
    if region is not None:
        cand &= region
    labels, n = ndimage.label(cand)
    if n == 0:
        return 0
    sizes = ndimage.sum(cand, labels, index=np.arange(1, n + 1))
    return int((sizes >= min_area_px).sum())
