"""Pre-specified erasure definitions (paper Section 3.7).

Thresholds z_det, z_miss and the interval level are read from the frozen
pre-registration config, never chosen after the data are seen.
"""
from __future__ import annotations

import numpy as np


def is_erased(z_ref: float, z_recon: float, z_det: float, z_miss: float) -> bool:
    """Detectable in the fully sampled reference, not in the reconstruction."""
    return bool(z_ref >= z_det and z_recon < z_miss)


def within_central_interval(value: float, samples: np.ndarray, level: float = 0.95) -> bool:
    lo, hi = np.quantile(samples, [(1 - level) / 2, 1 - (1 - level) / 2])
    return bool(lo <= value <= hi)


def is_silently_erased(z_ref: float, z_recon: float, psnr_present: float,
                       psnr_absent_samples: np.ndarray, z_det: float, z_miss: float,
                       level: float = 0.95) -> bool:
    """Erased, and global PSNR sits inside the lesion-free central interval."""
    return is_erased(z_ref, z_recon, z_det, z_miss) and within_central_interval(
        psnr_present, psnr_absent_samples, level)
