"""Classical reconstructors used for tests, smoke runs and as reference points.

Every reconstructor is a callable (y_masked, mask, maps, seed) -> complex image,
the interface the audit runner expects from the learned models too.
"""
from __future__ import annotations

import numpy as np

from csa.nullspace.projectors import cg_solve
from csa.physics.operator import SenseOperator


class ZeroFilled:
    name = "zero_filled"

    def __call__(self, y_masked, mask, maps, seed=0):
        return SenseOperator(maps, mask).adjoint(y_masked)


class CGSense:
    """Tikhonov-regularised SENSE: (A^H A + lam I)^{-1} A^H y by conjugate gradients.

    Linear and data-consistent, so for a pair (Y, Y + A ell) the lesion transfer
    is exactly P_lam ell when the same lam is used: t_R = 1 by construction.
    """
    name = "cg_sense"

    def __init__(self, lam: float = 1e-3, tol: float = 1e-8, maxiter: int = 500):
        self.lam, self.tol, self.maxiter = lam, tol, maxiter

    def __call__(self, y_masked, mask, maps, seed=0):
        op = SenseOperator(maps, mask)
        b = op.adjoint(y_masked)
        x, _ = cg_solve(lambda v: op.normal(v) + self.lam * v, b, tol=self.tol, maxiter=self.maxiter)
        return x
