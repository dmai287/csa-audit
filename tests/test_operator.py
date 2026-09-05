"""Gate 2: the adjoint is the true adjoint (dot-product test)."""
import numpy as np

from csa.physics.operator import SenseOperator, fft2c, ifft2c


def test_fft_is_unitary():
    rng = np.random.default_rng(0)
    x = rng.standard_normal((32, 32)) + 1j * rng.standard_normal((32, 32))
    assert np.allclose(ifft2c(fft2c(x)), x, atol=1e-12)
    assert np.isclose(np.linalg.norm(fft2c(x)), np.linalg.norm(x))


def test_dot_product_identity(scene):
    op = scene["op"]
    rng = np.random.default_rng(1)
    n = scene["n"]
    x = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    y = rng.standard_normal((op.n_coils, n, n)) + 1j * rng.standard_normal((op.n_coils, n, n))
    lhs = np.vdot(op.forward(x), y)
    rhs = np.vdot(x, op.adjoint(y))
    assert abs(lhs - rhs) / abs(lhs) < 1e-10


def test_mask_shapes(scene):
    op = SenseOperator(scene["maps"], scene["mask"])
    assert op.mask.shape == scene["maps"].shape[-2:]
    assert op.forward(scene["X"]).shape == (op.n_coils,) + scene["maps"].shape[-2:]
