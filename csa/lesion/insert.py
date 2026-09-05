"""Algorithm 1: additive lesion insertion carried through the forward operator.

The lesion is an image-domain perturbation ell with compact support, a
magnitude set by a contrast relative to the local tissue signal, and the phase
of the background so that it adds coherently. Because A is linear,
Y_cf = Y + A ell reuses the measured noise, mask and coil maps exactly.
"""
from __future__ import annotations

import numpy as np

from csa.physics.operator import SenseOperator, fft2c


def sphere_radius_mm(volume_mm3: float) -> float:
    """Radius of a sphere with the given volume."""
    return (3.0 * volume_mm3 / (4.0 * np.pi)) ** (1.0 / 3.0)


def slice_fill_fraction(radius_mm: float, thickness_mm: float) -> float:
    """Partial-volume factor for a sphere centred in a slice of given thickness.

    The in-plane disc of the sphere's equator is scaled by the ratio of the
    sphere's volume to the cylinder (disc area x slice thickness); capped at 1.
    """
    return float(min(1.0, 4.0 * radius_mm / (3.0 * thickness_mm)))


def soft_disc(shape, center, radius_px: float, taper_px: float = 1.0) -> np.ndarray:
    """Disc of radius `radius_px` at `center` (row, col) with a linear edge taper."""
    rr, cc = np.mgrid[0:shape[0], 0:shape[1]]
    r = np.hypot(rr - center[0], cc - center[1])
    if taper_px <= 0:
        return (r <= radius_px).astype(float)
    return np.clip((radius_px + taper_px / 2.0 - r) / taper_px, 0.0, 1.0)


def lesion_perturbation(X: np.ndarray, center, radius_px: float, contrast: float,
                        fill: float = 1.0, taper_px: float = 1.0) -> np.ndarray:
    """Complex perturbation ell: magnitude contrast*fill*mean|X|_disc, phase angle(X).

    contrast > 0 gives a hyperintense lesion, contrast < 0 a hypointense one.
    """
    m = soft_disc(X.shape, center, radius_px, taper_px)
    core = m > 0.5
    a = np.abs(X)
    ref = float(a[core].mean()) if core.any() else float(a[tuple(center)])
    if ref <= 0.0:
        # local tissue has no signal (e.g. CSF in a phantom): reference the surrounding tissue
        rr, cc = np.mgrid[0:X.shape[0], 0:X.shape[1]]
        r = np.hypot(rr - center[0], cc - center[1])
        ring = (r > radius_px + 1) & (r <= radius_px + 6) & (a > 0)
        ref = float(a[ring].mean()) if ring.any() else float(a[a > 0].mean())
    if ref <= 0.0:
        raise ValueError("cannot set lesion contrast: no tissue signal near the site")
    magnitude = contrast * fill * ref * m
    phase = np.angle(X)
    return magnitude * np.exp(1j * phase)


def counterfactual_kspace(op: SenseOperator, y_masked: np.ndarray, ell: np.ndarray) -> np.ndarray:
    """Y_cf = Y + A ell (Eq. additive_cf)."""
    return y_masked + op.forward(ell)


def insert_fully_sampled(maps: np.ndarray, y_full: np.ndarray, ell: np.ndarray) -> np.ndarray:
    """Alternative path: add the lesion to the fully sampled coil k-space before masking."""
    return y_full + fft2c(maps * ell[None])
