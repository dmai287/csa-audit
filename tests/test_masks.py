import numpy as np

from csa.physics.masks import equispaced_mask, random_mask


def test_equispaced_fraction_close_to_one_over_R():
    for R, cf in [(4, 0.08), (6, 0.06), (8, 0.04)]:
        m = equispaced_mask(320, R, cf)
        assert abs(m.mean() - 1.0 / R) < 0.02


def test_random_mask_seeded_and_centre_sampled():
    a, b = random_mask(320, 8, 0.04, seed=3), random_mask(320, 8, 0.04, seed=3)
    assert np.array_equal(a, b)
    assert a[160 - 6:160 + 6].all()
