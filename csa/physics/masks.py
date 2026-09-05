"""Cartesian undersampling masks following the fastMRI conventions.

`equispaced_mask` reproduces fastMRI's EquispacedMaskFractionFunc: a fully
sampled centre of `center_fraction` lines plus equispaced lines whose spacing
is adjusted so the overall sampled fraction is 1/acceleration.
"""
from __future__ import annotations

import numpy as np


def _num_low(n_cols: int, center_fraction: float) -> int:
    return int(round(n_cols * center_fraction))


def center_mask(n_cols: int, center_fraction: float) -> np.ndarray:
    mask = np.zeros(n_cols, dtype=bool)
    k = _num_low(n_cols, center_fraction)
    pad = (n_cols - k + 1) // 2
    mask[pad:pad + k] = True
    return mask


def equispaced_mask(n_cols: int, acceleration: int, center_fraction: float,
                    offset: int = 0) -> np.ndarray:
    """(W,) boolean mask along the phase-encode axis."""
    k = _num_low(n_cols, center_fraction)
    mask = center_mask(n_cols, center_fraction)
    adjusted = (acceleration * (k - n_cols)) / (k * acceleration - n_cols)
    samples = np.arange(offset, n_cols - 1, adjusted)
    mask[np.around(samples).astype(int)] = True
    return mask


def random_mask(n_cols: int, acceleration: int, center_fraction: float,
                seed: int = 0) -> np.ndarray:
    """(W,) boolean mask: fully sampled centre plus uniformly random lines."""
    k = _num_low(n_cols, center_fraction)
    prob = (n_cols / acceleration - k) / (n_cols - k)
    rng = np.random.default_rng(seed)
    mask = rng.uniform(size=n_cols) < prob
    mask |= center_mask(n_cols, center_fraction)
    return mask


def default_center_fraction(acceleration: int) -> float:
    """fastMRI pairs R=4 with 0.08 and R=8 with 0.04; R=6 is interpolated."""
    table = {4: 0.08, 6: 0.06, 8: 0.04}
    return table.get(int(acceleration), 0.32 / acceleration)
