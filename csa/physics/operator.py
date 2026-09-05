"""The multi-coil forward operator A = M F S and its adjoint.

Centred, orthonormal FFTs are used throughout so that `adjoint` is the exact
adjoint of `forward` (tests/test_operator.py checks the dot-product identity).
"""
from __future__ import annotations

import numpy as np


def fft2c(x: np.ndarray) -> np.ndarray:
    """Centred orthonormal 2-D FFT over the last two axes."""
    return np.fft.fftshift(
        np.fft.fft2(np.fft.ifftshift(x, axes=(-2, -1)), axes=(-2, -1), norm="ortho"),
        axes=(-2, -1),
    )


def ifft2c(y: np.ndarray) -> np.ndarray:
    """Centred orthonormal 2-D inverse FFT over the last two axes."""
    return np.fft.fftshift(
        np.fft.ifft2(np.fft.ifftshift(y, axes=(-2, -1)), axes=(-2, -1), norm="ortho"),
        axes=(-2, -1),
    )


def sense_combine(coil_images: np.ndarray, maps: np.ndarray) -> np.ndarray:
    """Sensitivity-weighted combination sum_c conj(S_c) x_c of coil images (C, H, W)."""
    return (np.conj(maps) * coil_images).sum(axis=0)


class SenseOperator:
    """A x = M * F (S_c x) per coil; A^H y = sum_c conj(S_c) F^{-1}(M * y_c).

    Parameters
    ----------
    maps : (C, H, W) complex coil sensitivities.
    mask : (H, W) boolean, or (W,) boolean for Cartesian phase-encode lines
           along the last axis.
    """

    def __init__(self, maps: np.ndarray, mask: np.ndarray):
        self.maps = np.asarray(maps)
        if self.maps.ndim != 3:
            raise ValueError("maps must have shape (C, H, W)")
        mask = np.asarray(mask, dtype=bool)
        if mask.ndim == 1:
            mask = np.broadcast_to(mask[None, :], self.maps.shape[-2:]).copy()
        if mask.shape != self.maps.shape[-2:]:
            raise ValueError("mask must have shape (H, W) or (W,)")
        self.mask = mask

    @property
    def n_coils(self) -> int:
        return self.maps.shape[0]

    @property
    def image_shape(self):
        return self.maps.shape[-2:]

    def forward(self, x: np.ndarray) -> np.ndarray:
        return self.mask[None] * fft2c(self.maps * x[None])

    def adjoint(self, y: np.ndarray) -> np.ndarray:
        return (np.conj(self.maps) * ifft2c(self.mask[None] * y)).sum(axis=0)

    def normal(self, x: np.ndarray) -> np.ndarray:
        """A^H A x."""
        return self.adjoint(self.forward(x))

    def relative_residual(self, x: np.ndarray, y: np.ndarray) -> float:
        """||A x - y||_2 / ||y||_2, the Phase 1 statistic."""
        return float(np.linalg.norm(self.forward(x) - y) / np.linalg.norm(y))

    def undersample(self, y_full: np.ndarray) -> np.ndarray:
        """Apply the mask to fully sampled multi-coil k-space."""
        return self.mask[None] * y_full
