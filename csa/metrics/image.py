"""Global and ROI image measures. PSNR and SSIM follow fastmri.evaluate conventions
(maxval = max of the reference, SSIM with a 7-pixel window and data_range = maxval)."""
from __future__ import annotations

import numpy as np
from skimage.metrics import structural_similarity


def _mag(x):
    return np.abs(np.asarray(x)).astype(float)


def psnr(gt: np.ndarray, pred: np.ndarray, maxval: float = None) -> float:
    gt, pred = _mag(gt), _mag(pred)
    maxval = float(gt.max()) if maxval is None else maxval
    mse = float(np.mean((gt - pred) ** 2))
    if mse == 0.0:
        return float("inf")
    return float(20 * np.log10(maxval) - 10 * np.log10(mse))


def ssim(gt: np.ndarray, pred: np.ndarray, maxval: float = None) -> float:
    gt, pred = _mag(gt), _mag(pred)
    maxval = float(gt.max()) if maxval is None else maxval
    return float(structural_similarity(gt, pred, data_range=maxval))


def roi_cnr(img: np.ndarray, roi: np.ndarray, background: np.ndarray) -> float:
    """(mean over roi - mean over background) / sd over background, on magnitude."""
    a = _mag(img)
    return float((a[roi].mean() - a[background].mean()) / a[background].std())


def roi_ssim(gt: np.ndarray, pred: np.ndarray, center, size: int, maxval: float = None) -> float:
    """SSIM restricted to a size x size window around (row, col)."""
    gt, pred = _mag(gt), _mag(pred)
    maxval = float(gt.max()) if maxval is None else maxval
    h = size // 2
    r0, c0 = max(int(center[0]) - h, 0), max(int(center[1]) - h, 0)
    g = gt[r0:r0 + size, c0:c0 + size]
    p = pred[r0:r0 + size, c0:c0 + size]
    win = min(7, (min(g.shape) // 2) * 2 - 1)
    return float(structural_similarity(g, p, data_range=maxval, win_size=win))


def center_crop(x: np.ndarray, size: int = 320) -> np.ndarray:
    """fastMRI's evaluation crop: the central size x size region (brain challenge convention).
    Global PSNR and SSIM are reported on this crop when `size` is set."""
    x = np.asarray(x)
    H, W = x.shape[-2:]
    h, w = min(size, H), min(size, W)
    r0, c0 = (H - h) // 2, (W - w) // 2
    return x[..., r0:r0 + h, c0:c0 + w]


def annulus(shape, center, r_in: float, r_out: float) -> np.ndarray:
    """Boolean annulus used as the local background for ROI CNR."""
    rr, cc = np.mgrid[0:shape[0], 0:shape[1]]
    r = np.hypot(rr - center[0], cc - center[1])
    return (r >= r_in) & (r < r_out)
