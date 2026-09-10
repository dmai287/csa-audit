"""Adapter for the csgm-mri-langevin reconstructor (Jalal et al 2021, NeurIPS).

The released sampler (annealed Langevin dynamics with an NCSNv2 score model
trained on fastMRI brain MVUE images) expects data prepared as in its own
dataloader: the readout field of view halved (fastMRI's readout oversampling
removed) and the phase-encode width resized to the model's image size (384
for the T2 config, 320 for FLAIR and T1), coil maps processed identically,
the masked multi-coil k-space as measurement, and the mask along the
phase-encode axis. This adapter performs exactly that preprocessing with the
same centred orthonormal transforms, calls the released `LangevinOptimizer`
unchanged, and maps the result back onto the audit's native grid: the
phase-encode padding is removed in k-space, and the readout margin outside
the model's field of view is filled with the zero-filled reconstruction
(measured data; no anatomy lives there) so that global scores stay
comparable across models.

Seeds: the sampler draws its Langevin noise from torch's global generator;
`seed` is applied immediately before sampling so both arms of a pair share
the same noise sequence (common random numbers), and the replicate seeds are
applied the same way on the replicate subsample.

The sampler is GPU-only as released (it casts to `torch.cuda.FloatTensor`).
The preprocessing and grid mapping are pure NumPy and are unit-tested on CPU
(`self_test`); the sampling call is not exercised without a GPU.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from csa.physics.operator import SenseOperator, fft2c, ifft2c

CONFIG_BY_SEQUENCE = {"AXT2": "brain_T2.yaml", "AXFLAIR": "brain_FLAIR.yaml", "AXT1": "brain_T1.yaml",
                      "AXT1POST": "brain_T1.yaml"}


def _resize_center(x: np.ndarray, size: int, axis: int) -> np.ndarray:
    """Centre crop or zero-pad `x` along `axis` to `size` (sigpy.resize semantics)."""
    n = x.shape[axis]
    if size == n:
        return x
    out_shape = list(x.shape); out_shape[axis] = size
    out = np.zeros(out_shape, dtype=x.dtype)
    if size < n:
        i = int(np.ceil((n - size) / 2.0))
        sl = [slice(None)] * x.ndim; sl[axis] = slice(i, i + size)
        return x[tuple(sl)].copy()
    i = int(np.ceil((size - n) / 2.0))
    sl = [slice(None)] * x.ndim; sl[axis] = slice(i, i + n)
    out[tuple(sl)] = x
    return out


def _fft_axis(x, axis):
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x, axes=axis), axis=axis, norm="ortho"), axes=axis)


def _ifft_axis(x, axis):
    return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(x, axes=axis), axis=axis, norm="ortho"), axes=axis)


def preprocess(kspace: np.ndarray, maps: np.ndarray, mask: np.ndarray, image_size):
    """Replicates csgm-mri-langevin's dataloader: (C, H, W) k-space and maps -> model grid.

    Returns k-space (C, h, w), maps (C, h, w), mask (w,), and the readout crop
    offset needed to map back. h = image_size[0] along readout (axis -2),
    w = image_size[1] along phase-encode (axis -1).
    """
    h, w = int(image_size[0]), int(image_size[1])
    # phase-encode: crop/pad in k-space
    ksp = _resize_center(kspace, w, axis=-1)
    # readout: halve the field of view in image space, back to k-space
    ksp = _fft_axis(_resize_center(_ifft_axis(ksp, -2), h, axis=-2), -2)
    # maps: same operations through k-space
    m = _resize_center(fft2c(maps), w, axis=-1)
    m = _fft_axis(_resize_center(_ifft_axis(m, -2), h, axis=-2), -2)
    m = ifft2c(m)
    mk = _resize_center(np.asarray(mask, dtype=float), w, axis=-1)
    ro_offset = int(np.ceil((kspace.shape[-2] - h) / 2.0))
    return ksp.astype(np.complex64), m.astype(np.complex64), mk, ro_offset


def mvue(kspace: np.ndarray, maps: np.ndarray) -> np.ndarray:
    """csgm-mri-langevin's MVUE estimate: sum(ifft(k) conj(S)) / sqrt(sum |S|^2)."""
    num = (ifft2c(kspace) * np.conj(maps)).sum(0)
    den = np.sqrt((np.abs(maps) ** 2).sum(0)); den[den == 0] = 1.0
    return num / den


def embed_back(x_model: np.ndarray, native_shape, image_size, ro_offset: int, native_pe: int) -> np.ndarray:
    """Model-grid image (h, w) -> native-grid image (H, W) on the model's readout FOV.

    Phase-encode: transform to k-space, crop/pad to the native width, back to
    image space. Readout: place into the native grid at `ro_offset`.
    """
    H, W = native_shape
    kx = _fft_axis(x_model, -1)
    kx = _resize_center(kx, W, axis=-1)
    x_pe = _ifft_axis(kx, -1)
    out = np.zeros((H, W), dtype=np.complex128)
    h = x_pe.shape[0]
    out[ro_offset:ro_offset + h, :] = x_pe
    return out


class DiffusionReconstructor:
    name = "diffusion_langevin"

    def __init__(self, repo_dir: str, checkpoint: str, sequence: str = "AXT2", device: str = "cuda",
                 config_overrides: dict = None):
        """repo_dir: a clone of https://github.com/utcsilab/csgm-mri-langevin; checkpoint: the
        NCSNv2 'mri-mvue' checkpoint_100000.pth; sequence selects the released config."""
        import yaml
        try:
            import torch
        except ImportError as exc:  # pragma: no cover
            raise ImportError("the diffusion adapter needs torch: pip install -r requirements-gpu.txt") from exc
        if not torch.cuda.is_available():
            raise RuntimeError("csgm-mri-langevin's sampler is GPU-only as released; no CUDA device found")
        self.torch = torch
        self.repo_dir = repo_dir
        if repo_dir not in sys.path:
            sys.path.insert(0, repo_dir)
        from main import LangevinOptimizer  # noqa: E402  (csgm-mri-langevin)
        cfg_path = os.path.join(repo_dir, "configs", "file", CONFIG_BY_SEQUENCE[sequence])
        cfg = yaml.safe_load(open(cfg_path))
        cfg["gen_ckpt"] = checkpoint; cfg["device"] = device; cfg["save_images"] = False; cfg["save_latent"] = False
        for k, v in (config_overrides or {}).items():
            cfg[k] = v
        self.cfg = cfg
        self.image_size = tuple(int(v) for v in cfg["image_size"])
        import logging
        self.sampler = LangevinOptimizer(cfg, logging.getLogger("csgm"), project_dir=repo_dir)

    def __call__(self, y_masked: np.ndarray, mask: np.ndarray, maps: np.ndarray, seed: int = 0) -> np.ndarray:
        torch = self.torch
        ksp, m, mk, ro_off = preprocess(y_masked, maps, mask, self.image_size)
        est = mvue(ksp.astype(np.complex128), m.astype(np.complex128))
        dev = self.cfg["device"]
        ref = torch.from_numpy(ksp)[None].to(dev).type(torch.complex128)
        maps_t = torch.from_numpy(m)[None].to(dev)
        mask_t = torch.from_numpy(mk)[None].to(dev)
        mvue_t = torch.from_numpy(est.astype(np.complex64))[None].to(dev)
        torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
        out = self.sampler.sample((ref, mvue_t, maps_t, mask_t))[0]["mvue"]  # (1, 2, h, w) normalised to the mvue scale
        o = out[0].detach().cpu().numpy()
        x_model = (o[0] + 1j * o[1]).astype(np.complex128)
        # back to the native grid; fill the readout margin with the zero-filled reconstruction
        native = SenseOperator(maps, mask).adjoint(y_masked)
        x = embed_back(x_model, native.shape, self.image_size, ro_off, y_masked.shape[-1])
        h = self.image_size[0]
        full = native.copy(); full[ro_off:ro_off + h, :] = x[ro_off:ro_off + h, :]
        return full


def self_test(kspace: np.ndarray, maps: np.ndarray, mask: np.ndarray, image_size=(320, 320)) -> dict:
    """CPU check of the preprocessing and the grid mapping, without the model."""
    ksp, m, mk, off = preprocess(kspace, maps, mask, image_size)
    est = mvue(ksp.astype(np.complex128), m.astype(np.complex128))
    # the native MVUE restricted to the readout FOV should match the model-grid MVUE mapped back
    native_mvue = mvue(kspace, maps)
    back = embed_back(est, native_mvue.shape, image_size, off, kspace.shape[-1])
    h = image_size[0]
    a, b = np.abs(back[off:off + h]), np.abs(native_mvue[off:off + h])
    corr = float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
    return {"model_grid": ksp.shape, "mask_width": int(mk.size), "sampled_fraction": float(mk.mean()),
            "readout_offset": off, "mvue_corr_native_vs_roundtrip": corr,
            "maps_norm_max": float(np.sqrt((np.abs(m) ** 2).sum(0)).max())}
