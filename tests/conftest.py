import numpy as np
import pytest

from csa.physics.masks import equispaced_mask
from csa.physics.operator import SenseOperator
from csa.synthetic import coil_maps, multicoil_kspace, phantom, tissue_site


@pytest.fixture(scope="session")
def scene():
    n, C = 96, 6
    X = phantom(n)
    maps = coil_maps(n, C)
    mask = equispaced_mask(n, acceleration=4, center_fraction=0.08)
    op = SenseOperator(maps, mask)
    y_full = multicoil_kspace(X, maps, noise_sigma=0.0)
    site = tissue_site(X, seed=0)
    assert np.abs(X)[site] > 0
    return {"n": n, "X": X, "maps": maps, "mask": mask, "op": op, "y_full": y_full, "site": site}
