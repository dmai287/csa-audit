"""Two-channel (real, imaginary) image-domain U-Net baseline.

The public fastMRI U-Net baseline is magnitude-only; the audit needs complex
images, so this baseline is trained with in_chans=out_chans=2 on the standard
fastMRI brain split (see scripts/train_unet_complex.py, to be written). Only
the inference adapter lives here.
"""
from __future__ import annotations

import numpy as np

from csa.physics.operator import SenseOperator


class ComplexUNetReconstructor:
    name = "unet_complex"

    def __init__(self, checkpoint: str, device: str = "cpu", chans: int = 32, num_pool_layers: int = 4):
        try:
            import torch
            from fastmri.models import Unet
        except ImportError as exc:  # pragma: no cover
            raise ImportError("U-Net needs torch and fastmri: pip install 'csa-audit[recon]'") from exc
        self.torch = torch
        self.device = torch.device(device)
        self.model = Unet(in_chans=2, out_chans=2, chans=chans, num_pool_layers=num_pool_layers)
        state = torch.load(checkpoint, map_location="cpu")
        self.model.load_state_dict(state.get("state_dict", state))
        self.model.eval().to(self.device)

    def __call__(self, y_masked: np.ndarray, mask: np.ndarray, maps: np.ndarray, seed: int = 0) -> np.ndarray:
        torch = self.torch
        torch.manual_seed(seed)
        zf = SenseOperator(maps, mask).adjoint(y_masked)
        scale = np.abs(zf).max() or 1.0
        x = torch.from_numpy(np.stack([zf.real, zf.imag], 0) / scale).float()[None].to(self.device)
        with torch.no_grad():
            out = self.model(x)[0].cpu().numpy() * scale
        return (out[0] + 1j * out[1]).astype(np.complex128)
