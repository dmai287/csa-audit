"""Both arms of a pair must be scored on a common intensity scale.

PSNR and SSIM take their scale from the reference maximum. The counterfactual
reference X + ell has a larger maximum whenever the lesion is bright, so
scoring each arm against its own reference shifts the score by
20 log10(max_cf / max_f) with no change in the error at all. These tests pin
the convention that the audit uses the lesion-absent reference's maximum for
both arms.
"""
import numpy as np

from csa.metrics.image import center_crop, psnr, ssim


def _scene(seed=0, n=96):
    rng = np.random.default_rng(seed)
    X = np.abs(rng.normal(1.0, 0.05, (n, n)))
    ell = np.zeros_like(X)
    ell[40:44, 40:44] = 1.5 * X.max()          # a bright lesion that raises the maximum
    err = rng.normal(0, 0.01, (n, n))           # identical error in both arms
    return X, ell, err


def test_own_maxval_shifts_psnr_without_any_error_change():
    X, ell, err = _scene()
    X_cf = X + ell
    p_own = psnr(X_cf, X_cf + err) - psnr(X, X + err)
    assert abs(p_own) > 0.5, "the artefact should be large for a bright lesion"
    shift = 20 * np.log10(X_cf.max() / X.max())
    assert np.isclose(p_own, shift, atol=1e-6)


def test_common_maxval_removes_the_artefact():
    X, ell, err = _scene()
    X_cf = X + ell
    m = float(X.max())
    d = psnr(X_cf, X_cf + err, maxval=m) - psnr(X, X + err, maxval=m)
    assert abs(d) < 1e-9, "with a common scale, an identical error must give an identical score"
    ds = ssim(X_cf, X_cf + err, maxval=m) - ssim(X, X + err, maxval=m)
    assert abs(ds) < 0.05


def test_center_crop_shape_and_centering():
    x = np.arange(640 * 320, dtype=float).reshape(640, 320)
    c = center_crop(x, 320)
    assert c.shape == (320, 320)
    assert np.allclose(c, x[160:480, :320])
