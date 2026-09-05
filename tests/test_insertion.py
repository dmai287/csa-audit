"""Gate 4: both insertion paths agree to machine precision; noise is shared exactly."""
import numpy as np

from csa.lesion.insert import (counterfactual_kspace, insert_fully_sampled, lesion_perturbation,
                               slice_fill_fraction, sphere_radius_mm)
from csa.physics.operator import SenseOperator, sense_combine, ifft2c
from csa.synthetic import multicoil_kspace


def test_two_insertion_paths_agree(scene):
    op, X, maps = scene["op"], scene["X"], scene["maps"]
    y_full = multicoil_kspace(X, maps, noise_sigma=0.05, seed=7)
    ell = lesion_perturbation(X, center=scene["site"], radius_px=2.5, contrast=0.4)
    y = op.undersample(y_full)
    path_a = counterfactual_kspace(op, y, ell)
    path_b = op.undersample(insert_fully_sampled(maps, y_full, ell))
    assert np.allclose(path_a, path_b, atol=1e-12)
    # the difference between the arms is exactly the lesion's projection: noise cancels
    assert np.allclose(path_a - y, op.forward(ell), atol=1e-12)


def test_fully_sampled_roundtrip_recovers_lesion(scene):
    X, maps = scene["X"], scene["maps"]
    full = SenseOperator(maps, np.ones(X.shape, dtype=bool))
    ell = lesion_perturbation(X, center=scene["site"], radius_px=3.0, contrast=0.5)
    y_cf = insert_fully_sampled(maps, multicoil_kspace(X, maps), ell)
    rec = sense_combine(ifft2c(y_cf), maps)
    assert np.allclose(rec, X + ell, atol=1e-10)


def test_lesion_phase_matches_background(scene):
    X = scene["X"]
    ell = lesion_perturbation(X, center=scene["site"], radius_px=3.0, contrast=0.5, taper_px=0.0)
    core = np.abs(ell) > 0
    assert np.allclose(np.angle(ell[core]), np.angle(X[core]))


def test_geometry_helpers():
    assert np.isclose(sphere_radius_mm(4.0 / 3.0 * np.pi), 1.0)
    assert slice_fill_fraction(1.0, 5.0) == 4.0 / 15.0
    assert slice_fill_fraction(10.0, 5.0) == 1.0
