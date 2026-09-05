"""Adapter for the fastMRI End-to-End Variational Network (Sriram et al 2020).

Requires torch and the `fastmri` package plus a checkpoint (the public brain
leaderboard weights). Import is lazy so the NumPy core works without torch.
"""
from __future__ import annotations

import numpy as np


class VarNetReconstructor:
    name = "e2e_varnet"

    def __init__(self, checkpoint: str, device: str = "cpu", num_cascades: int = 12,
                 sens_chans: int = 8, sens_pools: int = 4, chans: int = 18, pools: int = 4):
        try:
            import torch
            from fastmri.models import VarNet
        except ImportError as exc:  # pragma: no cover
            raise ImportError("VarNet needs torch and fastmri: pip install 'csa-audit[recon]'") from exc
        self.torch = torch
        self.device = torch.device(device)
        self.model = VarNet(num_cascades=num_cascades, sens_chans=sens_chans, sens_pools=sens_pools,
                            chans=chans, pools=pools)
        state = torch.load(checkpoint, map_location="cpu")
        state = state.get("state_dict", state)
        state = {k.replace("varnet.", "", 1) if k.startswith("varnet.") else k: v for k, v in state.items()}
        self.model.load_state_dict(state)
        self.model.eval().to(self.device)

    def __call__(self, y_masked: np.ndarray, mask: np.ndarray, maps: np.ndarray, seed: int = 0) -> np.ndarray:
        """Returns the VarNet output. NOTE: the public model returns a magnitude
        image (RSS of the coil estimates); complex output requires exposing the
        cascade output, which the audit needs for the null-space statistics.
        See docs/MODELS.md before using this adapter in an experiment."""
        torch = self.torch
        torch.manual_seed(seed)
        y = torch.from_numpy(np.stack([y_masked.real, y_masked.imag], -1)).float()[None].to(self.device)
        m = torch.from_numpy(np.asarray(mask, dtype=bool))
        m = m.reshape(1, 1, 1, -1, 1) if m.ndim == 1 else m[None, None, :, :, None]
        with torch.no_grad():
            out = self.model(y, m.to(self.device))
        return out[0].cpu().numpy().astype(np.complex128)
