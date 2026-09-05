"""Synthetic phantom, coil maps and noise for tests and smoke runs.

Nothing here touches fastMRI data. The phantom is the Shepp-Logan image with a
smooth phase; coil maps are smooth complex bumps normalised so that
sum_c |S_c|^2 = 1, which mirrors what ESPIRiT returns inside the support.
"""
from __future__ import annotations

import numpy as np
from skimage.data import shepp_logan_phantom
from skimage.transform import resize

from csa.physics.operator import fft2c


def phantom(n: int = 128, seed: int = 0) -> np.ndarray:
    """Complex-valued (n, n) phantom in [0, 1] magnitude with a smooth phase."""
    img = resize(shepp_logan_phantom(), (n, n), anti_aliasing=True)
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n] / n
    phase = 0.8 * np.sin(2 * np.pi * (0.7 * xx + 0.3 * yy) + rng.uniform(0, 2 * np.pi))
    return (img * np.exp(1j * phase)).astype(np.complex128)


def coil_maps(n: int = 128, n_coils: int = 8, seed: int = 0) -> np.ndarray:
    """Smooth complex coil sensitivities, shape (n_coils, n, n), sum |S|^2 = 1."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n] / (n - 1)
    maps = np.zeros((n_coils, n, n), dtype=np.complex128)
    for c in range(n_coils):
        ang = 2 * np.pi * c / n_coils
        cx, cy = 0.5 + 0.6 * np.cos(ang), 0.5 + 0.6 * np.sin(ang)
        r2 = (xx - cx) ** 2 + (yy - cy) ** 2
        mag = np.exp(-r2 / (2 * 0.35 ** 2))
        ph = 2 * np.pi * (0.3 * (xx - cx) + 0.2 * (yy - cy)) + rng.uniform(0, 2 * np.pi)
        maps[c] = mag * np.exp(1j * ph)
    norm = np.sqrt((np.abs(maps) ** 2).sum(axis=0))
    return maps / norm


def multicoil_kspace(x: np.ndarray, maps: np.ndarray, noise_sigma: float = 0.0,
                     seed: int = 0) -> np.ndarray:
    """Fully sampled multi-coil k-space of image x with optional complex Gaussian noise."""
    y = fft2c(maps * x[None])
    if noise_sigma > 0:
        rng = np.random.default_rng(seed)
        noise = rng.standard_normal(y.shape) + 1j * rng.standard_normal(y.shape)
        y = y + noise_sigma * noise / np.sqrt(2.0)
    return y


def tissue_site(x: np.ndarray, seed: int = 0, lo: float = 0.3, hi: float = 0.9):
    """A (row, col) inside tissue: magnitude between lo and hi times max|x|."""
    a = np.abs(x)
    ok = np.argwhere((a > lo * a.max()) & (a < hi * a.max()))
    rng = np.random.default_rng(seed)
    r, c = ok[rng.integers(len(ok))]
    return int(r), int(c)


def support_mask(x: np.ndarray, frac: float = 0.05) -> np.ndarray:
    """Boolean mask of pixels whose magnitude exceeds frac * max|x|."""
    a = np.abs(x)
    return a > frac * a.max()
