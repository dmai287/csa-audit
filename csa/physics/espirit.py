"""Coil-sensitivity estimation. Real data: ESPIRiT via sigpy (optional dependency)."""
from __future__ import annotations

import numpy as np


def espirit_maps(kspace_full: np.ndarray, calib_width: int = 24, thresh: float = 0.02,
                 kernel_width: int = 6, crop: float = 0.95, max_iter: int = 100,
                 device: int = -1) -> np.ndarray:
    """ESPIRiT sensitivity maps (C, H, W) from fully sampled multi-coil k-space.

    Requires `sigpy`. Returns the first eigenvector map, as used to build the
    reference image and the counterfactual arm (paper Section 3.2).
    """
    try:
        import sigpy as sp
        import sigpy.mri as mr
    except ImportError as exc:  # pragma: no cover - exercised only without sigpy
        raise ImportError("ESPIRiT needs sigpy: pip install 'csa-audit[recon]'") from exc
    dev = sp.Device(device)
    ksp = sp.to_device(np.asarray(kspace_full, dtype=np.complex64), dev)
    app = mr.app.EspiritCalib(ksp, calib_width=calib_width, thresh=thresh,
                              kernel_width=kernel_width, crop=crop, max_iter=max_iter,
                              device=dev, show_pbar=False)
    return np.asarray(sp.to_device(app.run(), sp.cpu_device))
