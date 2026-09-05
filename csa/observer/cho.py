"""Channelised Hotelling observer with Laguerre-Gauss channels (SKE/BKS task).

Template w = K_v^{-1} (mean v1 - mean v0) is estimated on a training partition;
d', AUC and the per-lesion z-score are computed on a disjoint test partition.
"""
from __future__ import annotations

import numpy as np
from scipy.special import eval_laguerre
from scipy.stats import mannwhitneyu


def laguerre_gauss_channels(size: int, n_channels: int, a_px: float) -> np.ndarray:
    """(size*size, n_channels) matrix of rotationally symmetric LG channels.

    LG_j(r) = sqrt(2)/a * exp(-pi r^2/a^2) * L_j(2 pi r^2 / a^2), j = 0..n-1.
    """
    c = (size - 1) / 2.0
    yy, xx = np.mgrid[0:size, 0:size]
    r2 = (yy - c) ** 2 + (xx - c) ** 2
    u = 2.0 * np.pi * r2 / a_px ** 2
    U = np.zeros((size * size, n_channels))
    for j in range(n_channels):
        ch = (np.sqrt(2.0) / a_px) * np.exp(-u / 2.0) * eval_laguerre(j, u)
        U[:, j] = ch.ravel()
    return U


def extract_roi(img: np.ndarray, center, size: int) -> np.ndarray:
    """Magnitude patch of `size` x `size` centred at (row, col), zero-padded at borders."""
    a = np.abs(img)
    h = size // 2
    pad = np.pad(a, h, mode="constant")
    r, c = int(center[0]) + h, int(center[1]) + h
    return pad[r - h:r - h + size, c - h:c - h + size].astype(float)


class CHO:
    """Channelised Hotelling observer."""

    def __init__(self, channels: np.ndarray, ridge: float = 1e-8):
        self.U = channels
        self.ridge = ridge
        self.w = None

    def channel_outputs(self, rois: np.ndarray) -> np.ndarray:
        rois = np.asarray(rois, dtype=float).reshape(len(rois), -1)
        return rois @ self.U

    def fit(self, rois_present: np.ndarray, rois_absent: np.ndarray) -> "CHO":
        v1, v0 = self.channel_outputs(rois_present), self.channel_outputs(rois_absent)
        K = 0.5 * (np.cov(v1, rowvar=False) + np.cov(v0, rowvar=False))
        K = K + self.ridge * np.trace(K) / K.shape[0] * np.eye(K.shape[0])
        self.w = np.linalg.solve(K, v1.mean(0) - v0.mean(0))
        return self

    def decide(self, rois: np.ndarray) -> np.ndarray:
        if self.w is None:
            raise RuntimeError("call fit() first")
        return self.channel_outputs(rois) @ self.w


def dprime(t1: np.ndarray, t0: np.ndarray) -> float:
    """SNR of the decision variable between the two classes."""
    return float((np.mean(t1) - np.mean(t0)) / np.sqrt(0.5 * (np.var(t1, ddof=1) + np.var(t0, ddof=1))))


def auc_mann_whitney(t1: np.ndarray, t0: np.ndarray) -> float:
    """Non-parametric AUC = U / (n1 n0)."""
    u = mannwhitneyu(t1, t0, alternative="two-sided").statistic
    return float(u / (len(t1) * len(t0)))


def per_lesion_z(t1_i: float, t0: np.ndarray) -> float:
    """Standardised position of one lesion-present statistic in the lesion-absent distribution."""
    return float((t1_i - np.mean(t0)) / np.std(t0, ddof=1))


def ideal_dprime_white(signal: np.ndarray, sigma: float) -> float:
    """Ideal linear observer d' for a known signal in white Gaussian noise: ||s|| / sigma."""
    return float(np.linalg.norm(signal) / sigma)
