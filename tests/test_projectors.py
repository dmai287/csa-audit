"""Gate 5: filtered projectors match closed forms and the transfer identities hold."""
import numpy as np

from csa.lesion.insert import lesion_perturbation
from csa.nullspace.projectors import (cg_solve, measurement_fraction, measurement_gain,
                                      null_filter, range_filter, real_inner, transfer_statistics)
from csa.physics.masks import equispaced_mask
from csa.physics.operator import SenseOperator, ifft2c
from csa.recon.base import CGSense


def _single_coil_op(n=64, R=4):
    maps = np.ones((1, n, n), dtype=complex)
    return SenseOperator(maps, equispaced_mask(n, R, 0.08))


def test_single_coil_closed_form():
    """A^H A = F^H M F, so P_lam maps a sampled Fourier mode to x/(1+lam), an unsampled one to 0."""
    op = _single_coil_op()
    n = op.image_shape[0]
    lam = 0.3
    cols = np.where(op.mask[0])[0]
    uncols = np.where(~op.mask[0])[0]
    for col, expect in [(cols[len(cols) // 3], 1.0 / (1 + lam)), (uncols[len(uncols) // 2], 0.0)]:
        e = np.zeros((n, n), dtype=complex)
        e[n // 2 + 3, col] = 1.0
        x = ifft2c(e)
        p = range_filter(op, x, lam, tol=1e-12)
        assert np.allclose(p, expect * x, atol=1e-9)


def test_projector_is_hermitian(scene):
    op = scene["op"]
    rng = np.random.default_rng(5)
    n = scene["n"]
    x = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    y = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    lam = 1e-2
    lhs = np.vdot(range_filter(op, x, lam, tol=1e-12), y)
    rhs = np.vdot(x, range_filter(op, y, lam, tol=1e-12))
    assert abs(lhs - rhs) / abs(lhs) < 1e-8


def test_cg_solves_regularised_normal_equations(scene):
    op = scene["op"]
    rng = np.random.default_rng(2)
    n = scene["n"]
    b = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    lam = 1e-2
    z, info = cg_solve(lambda v: op.normal(v) + lam * v, b, tol=1e-12)
    assert info["converged"]
    assert np.linalg.norm(op.normal(z) + lam * z - b) / np.linalg.norm(b) < 1e-10


def test_transfer_identities(scene):
    op, X = scene["op"], scene["X"]
    ell = lesion_perturbation(X, center=scene["site"], radius_px=2.5, contrast=0.4)
    lam = 1e-3
    t_R, t_N = transfer_statistics(op, ell, ell, lam, tol=1e-12)
    assert np.isclose(t_R, 1.0, atol=1e-8) and np.isclose(t_N, 1.0, atol=1e-8)
    # a purely range-space transfer carries little null-space content
    _, t_N_range = transfer_statistics(op, range_filter(op, ell, lam, tol=1e-12), ell, lam, tol=1e-12)
    assert 0.0 <= t_N_range < 0.2
    mu = measurement_fraction(op, ell, lam, tol=1e-12)
    assert 0.0 < mu < 1.0
    assert 0.0 < measurement_gain(op, ell) < 1.0


def test_cg_sense_transfer_is_the_range_filter(scene):
    """A linear data-consistent reconstructor moves exactly P_lam ell between the arms."""
    op, X, maps, mask = scene["op"], scene["X"], scene["maps"], scene["mask"]
    lam = 1e-3
    ell = lesion_perturbation(X, center=scene["site"], radius_px=2.5, contrast=0.4)
    y = op.undersample(scene["y_full"])
    rec = CGSense(lam=lam, tol=1e-12)
    dX = rec(y + op.forward(ell), mask, maps) - rec(y, mask, maps)
    assert np.allclose(dX, range_filter(op, ell, lam, tol=1e-12), atol=1e-8)
    t_R, t_N = transfer_statistics(op, dX, ell, lam, tol=1e-12)
    # Tikhonov shrinkage of measured components costs O(lam) in t_R; hard DC would give exactly 1
    assert 0.98 < t_R <= 1.0 + 1e-9
    assert 0.0 <= t_N < 0.2


def test_residual_insensitivity_bound(scene):
    """Proposition 1 on a real operator: the residual change is bounded by ||A delta||."""
    op, X = scene["op"], scene["X"]
    y = op.undersample(scene["y_full"])
    rng = np.random.default_rng(9)
    n = scene["n"]
    Xh = X + 0.01 * (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
    delta = null_filter(op, lesion_perturbation(X, scene["site"], 2.5, 0.4), 1e-6, tol=1e-12)
    r0 = np.linalg.norm(op.forward(Xh) - y)
    r1 = np.linalg.norm(op.forward(Xh + delta) - y)
    assert abs(r1 - r0) <= np.linalg.norm(op.forward(delta)) + 1e-12
