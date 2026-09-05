"""Gate 6: the CHO never beats the ideal linear observer and scales with contrast."""
import numpy as np

from csa.observer.cho import (CHO, auc_mann_whitney, dprime, ideal_dprime_white,
                              laguerre_gauss_channels, per_lesion_z)


def _gaussian_blob(size, sigma, amp):
    c = (size - 1) / 2
    yy, xx = np.mgrid[0:size, 0:size]
    return amp * np.exp(-((yy - c) ** 2 + (xx - c) ** 2) / (2 * sigma ** 2))


def _run(amp, seed=0, size=24, sigma_n=1.0, n=600):
    rng = np.random.default_rng(seed)
    s = _gaussian_blob(size, 2.0, amp)
    absent = sigma_n * rng.standard_normal((n, size, size))
    present = absent + s  # same noise fields for a paired comparison
    U = laguerre_gauss_channels(size, 6, a_px=8.0)
    cho = CHO(U).fit(present[: n // 2].reshape(n // 2, -1), absent[: n // 2].reshape(n // 2, -1))
    t1 = cho.decide(present[n // 2:].reshape(n // 2, -1))
    t0 = cho.decide(absent[n // 2:].reshape(n // 2, -1))
    return dprime(t1, t0), ideal_dprime_white(s, sigma_n), t1, t0


def test_cho_bounded_by_ideal_and_linear_in_contrast():
    d1, ideal1, t1, t0 = _run(0.5)
    d2, ideal2, _, _ = _run(1.0)
    assert d1 <= 1.15 * ideal1
    assert d2 <= 1.15 * ideal2
    assert 0.7 * ideal1 < d1  # LG channels capture a Gaussian blob well
    assert 1.6 < d2 / d1 < 2.4
    auc = auc_mann_whitney(t1, t0)
    assert 0.5 < auc <= 1.0
    z = per_lesion_z(float(t1.mean()), t0)
    assert z > 0
