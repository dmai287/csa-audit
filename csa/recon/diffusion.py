"""Adapter for a score-based / Langevin MRI reconstructor with data consistency.

Target implementation: Jalal et al. (2021) csgm-mri-langevin, whose released
weights were trained on fastMRI brain multi-coil data. Common random numbers
are enforced by seeding the sampler identically for both arms of a pair.
"""
from __future__ import annotations

import numpy as np


class DiffusionReconstructor:
    name = "diffusion_langevin"

    def __init__(self, config_path: str, checkpoint: str, device: str = "cpu"):
        self.config_path, self.checkpoint, self.device = config_path, checkpoint, device
        raise NotImplementedError(
            "Wire this adapter to the csgm-mri-langevin sampler once its weights are "
            "downloaded; see docs/MODELS.md. Both arms must share the seed.")

    def __call__(self, y_masked: np.ndarray, mask: np.ndarray, maps: np.ndarray, seed: int = 0) -> np.ndarray:  # pragma: no cover
        raise NotImplementedError
