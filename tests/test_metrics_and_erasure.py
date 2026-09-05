import numpy as np

from csa.audit.erasure import is_erased, is_silently_erased, within_central_interval
from csa.metrics.image import annulus, psnr, roi_cnr, roi_ssim, ssim
from csa.stats.bootstrap import design_effect, pairs_per_cell, standardised_mean_difference


def test_psnr_ssim_identity(scene):
    X = np.abs(scene["X"])
    assert psnr(X, X) == np.inf or psnr(X, X) > 100
    assert np.isclose(ssim(X, X), 1.0)


def test_roi_measures(scene):
    X = np.abs(scene["X"])
    center = (48, 40)
    disc = np.hypot(*np.mgrid[0:96, 0:96] - np.array(center)[:, None, None]) <= 3
    ring = annulus(X.shape, center, 5, 12)
    bright = X.copy()
    bright[disc] += 0.5
    assert roi_cnr(bright, disc, ring) > roi_cnr(X, disc, ring)
    assert 0 < roi_ssim(X, bright, center, 16) < 1


def test_erasure_logic():
    assert is_erased(4.0, 1.0, z_det=3.0, z_miss=1.645)
    assert not is_erased(2.0, 1.0, z_det=3.0, z_miss=1.645)
    samples = np.linspace(30, 34, 101)
    assert within_central_interval(32.0, samples)
    assert not within_central_interval(35.0, samples)
    assert is_silently_erased(4.0, 1.0, 32.0, samples, 3.0, 1.645)
    assert not is_silently_erased(4.0, 1.0, 35.0, samples, 3.0, 1.645)


def test_sample_size_rule():
    assert np.isclose(design_effect(12, 0.1), 2.1)
    n = pairs_per_cell(pi0=0.05, half_width=0.02, m=12, rho=0.1)
    assert 900 < n < 1000
    assert np.isclose(standardised_mean_difference(np.array([1., 2., 3.]), np.array([1., 2., 3.])), 0.0)
