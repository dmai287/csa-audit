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

    def last_cascade_coil_images(self, y_masked: np.ndarray, mask: np.ndarray):
        """Run the cascades and return the network's multi-coil k-space and its coil images (C, H, W) complex."""
        torch, fastmri = self.torch, self.fastmri
        y = _to_torch_kspace(y_masked, torch).to(self.device)
        m = torch.from_numpy(np.asarray(mask, dtype=bool))
        m = (m.reshape(1, 1, 1, -1, 1) if m.ndim == 1 else m[None, None, :, :, None]).to(self.device)
        with torch.no_grad():
            sens = self.model.sens_net(y, m)               # (1, C, H, W, 2), the network's own maps
            kspace_pred = y.clone()
            for cascade in self.model.cascades:
                kspace_pred = cascade(kspace_pred, y, m, sens)
            coil_imgs = fastmri.ifft2c(kspace_pred)           # (1, C, H, W, 2)
        k = kspace_pred[0].cpu().numpy(); c = coil_imgs[0].cpu().numpy(); s_ = sens[0].cpu().numpy()
        return ((k[..., 0] + 1j * k[..., 1]).astype(np.complex128), (c[..., 0] + 1j * c[..., 1]).astype(np.complex128),
                (s_[..., 0] + 1j * s_[..., 1]).astype(np.complex128))

    def _complex_image(self, y_masked: np.ndarray, mask: np.ndarray, maps: np.ndarray = None):
        """SENSE-combine the network's coil images.

        With `maps` given (the audit's ESPIRiT maps, the default through the
        reconstructor interface) the combination is in the audit's phase
        convention, so the result is consistent with the audit's forward
        operator. Without them the network's own maps are used, which carry
        a different phase convention and only the magnitude is comparable.
        """
        _, coil_imgs, sens = self.last_cascade_coil_images(y_masked, mask)
        use = maps if (maps is not None and self.output == "complex_sense") else sens
        return (np.conj(use) * coil_imgs).sum(axis=0)

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
        return self._complex_image(y_masked, mask, maps)
