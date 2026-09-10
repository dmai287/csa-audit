"""Adapter for the fastMRI End-to-End Variational Network (Sriram et al 2020).

The public VarNet returns a root-sum-of-squares magnitude image. The audit
needs a complex image, so this adapter reproduces VarNet's forward pass up to
the final cascade output (multi-coil k-space), inverse-transforms it, and
combines the coils with the network's own estimated sensitivity maps. This
uses the released weights unchanged; only the read-out differs, and the
magnitude of the returned image equals the public output up to the
coil-combination rule (SENSE combination with the net's maps rather than RSS).

Weights: brain leaderboard checkpoint, see docs/MODELS.md. Import of torch and
fastmri is lazy so the NumPy core works without them.
"""
from __future__ import annotations

import numpy as np


def _to_torch_kspace(y: np.ndarray, torch):
    """(C, H, W) complex numpy -> (1, C, H, W, 2) float tensor."""
    return torch.from_numpy(np.stack([y.real, y.imag], -1).astype(np.float32))[None]


class VarNetReconstructor:
    name = "e2e_varnet"

    def __init__(self, checkpoint: str, device: str = "cuda", num_cascades: int = 12,
                 sens_chans: int = 8, sens_pools: int = 4, chans: int = 18, pools: int = 4,
                 output: str = "complex_sense"):
        try:
            import torch
            import fastmri
            from fastmri.models import VarNet
        except ImportError as exc:  # pragma: no cover
            raise ImportError("VarNet needs torch and fastmri: pip install -r requirements-gpu.txt") from exc
        self.torch, self.fastmri = torch, fastmri
        self.device = torch.device(device if (device != "cuda" or torch.cuda.is_available()) else "cpu")
        self.output = output
        self.model = VarNet(num_cascades=num_cascades, sens_chans=sens_chans, sens_pools=sens_pools, chans=chans, pools=pools)
        state = torch.load(checkpoint, map_location="cpu")
        state = state.get("state_dict", state)
        state = {k[len("varnet."):] if k.startswith("varnet.") else k: v for k, v in state.items()}
        self.model.load_state_dict(state)
        self.model.eval().to(self.device)

    def _complex_image(self, y_masked: np.ndarray, mask: np.ndarray):
        """Run the cascades and return the SENSE-combined complex image (H, W) as numpy."""
        torch, fastmri = self.torch, self.fastmri
        y = _to_torch_kspace(y_masked, torch).to(self.device)
        m = torch.from_numpy(np.asarray(mask, dtype=bool))
        m = (m.reshape(1, 1, 1, -1, 1) if m.ndim == 1 else m[None, None, :, :, None]).to(self.device)
        with torch.no_grad():
            sens = self.model.sens_net(y, m)               # (1, C, H, W, 2)
            kspace_pred = y.clone()
            for cascade in self.model.cascades:
                kspace_pred = cascade(kspace_pred, y, m, sens)
            coil_imgs = fastmri.ifft2c(kspace_pred)           # (1, C, H, W, 2)
            if self.output == "complex_sense":
                combined = fastmri.complex_mul(coil_imgs, fastmri.complex_conj(sens)).sum(dim=1)  # (1, H, W, 2)
            else:  # "rss_phase_of_first": magnitude by RSS, phase from the SENSE combination
                combined = fastmri.complex_mul(coil_imgs, fastmri.complex_conj(sens)).sum(dim=1)
                mag = fastmri.rss(fastmri.complex_abs(coil_imgs), dim=1)
                ph = torch.atan2(combined[..., 1], combined[..., 0])
                combined = torch.stack([mag * torch.cos(ph), mag * torch.sin(ph)], -1)
        c = combined[0].cpu().numpy()
        return (c[..., 0] + 1j * c[..., 1]).astype(np.complex128)

    def public_magnitude(self, y_masked: np.ndarray, mask: np.ndarray) -> np.ndarray:
        """The public VarNet output (RSS magnitude), for checking the adapter against the model as released."""
        torch = self.torch
        y = _to_torch_kspace(y_masked, torch).to(self.device)
        m = torch.from_numpy(np.asarray(mask, dtype=bool))
        m = (m.reshape(1, 1, 1, -1, 1) if m.ndim == 1 else m[None, None, :, :, None]).to(self.device)
        with torch.no_grad():
            out = self.model(y, m)
        return out[0].cpu().numpy()

    def __call__(self, y_masked: np.ndarray, mask: np.ndarray, maps: np.ndarray, seed: int = 0) -> np.ndarray:
        self.torch.manual_seed(seed)   # VarNet is deterministic; the seed is kept for interface uniformity
        return self._complex_image(y_masked, mask)
