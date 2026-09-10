# Models under audit

| Family | Implementation | Weights | Output | Status (2026-09-10) |
|---|---|---|---|---|
| E2E-VarNet | `csa/recon/varnet.py` over `fastmri.models.VarNet` | brain leaderboard checkpoint (`dl.fbaipublicfiles.com`, 114 MB) | complex: last-cascade k-space, inverse FFT, SENSE combination with the network's own maps; `public_magnitude()` returns the released RSS output for comparison | adapter written; CPU check on a real slice in `docs/VALIDATION.md` |
| Complex U-Net | `csa/recon/unet.py` over `fastmri.models.Unet(in_chans=2, out_chans=2)` | trained by `scripts/train_unet_complex.py` on train batch 0 (annotated slices and validation hosts excluded) | complex | training script written and CPU smoke-tested; ~1 GPU-day for the real model |
| Diffusion | `csa/recon/diffusion.py` over csgm-mri-langevin's `LangevinOptimizer` (Jalal et al 2021) | `checkpoint_100000.pth` (NCSNv2 'mri-mvue', 1.5 GB, Google Drive link in the repo README) | complex, on the released model grid (readout FOV halved; phase-encode 384 for the T2 config, 320 for FLAIR/T1), mapped back to the native grid with the zero-filled reconstruction filling the readout margin | preprocessing and grid mapping unit-tested on real slices (MVUE round trip 0.999); the sampling call is GPU-only as released and untested here |

Notes on the diffusion adapter: the released dataloader generates its
undersampling mask over the padded phase-encode width; this adapter pads the
native mask with zeros instead, so no line that was never acquired counts as
measured. Both arms of a pair get the same `seed` immediately before
sampling (common random numbers); replicate seeds run on the 10 % subsample
only (`configs/confirmatory_A1b.yaml`).

Rules: models are audited as released; no model is retrained on the lesion
bank; no lesion-bearing slice enters any training or calibration step; both
arms of a pair use the same seed; seeds are then varied across replicates.

Chung and Ye's score-MRI weights were trained on knee data and are therefore
not the first choice for a brain audit.

## What is still needed to run an experiment

1. **Environment.** This repo's NumPy core needs no GPU. The three
   reconstructors need `pip install -e ".[recon]"` (torch, fastmri, sigpy) and
   BART for the diffusion sampler's sensitivity-map step; none of these are
   installed in the machine's default Python 3.8 environment as of
   2026-09-08.
2. **Compute.** This machine has no CUDA GPU (an Intel UHD 630 and an AMD
   Radeon Pro 560X, 4 GB VRAM, Metal only). `fastmri`, VarNet training/
   inference and Langevin sampling assume CUDA; PyTorch's Metal (MPS) backend
   may run some of this but is unproven for these codebases and will be far
   slower. A CUDA machine or a cloud GPU instance is the practical requirement
   for Experiments 1 and 3, which reconstruct every pair with every model.
3. **The complex U-Net checkpoint**, trained from `train_batch_0`
   (`scripts/train_unet_complex.py` is not written yet).
4. **The VarNet adapter change** to expose the complex cascade output instead
   of the magnitude RSS image (`csa/recon/varnet.py`).
5. **The diffusion adapter**, still `NotImplementedError` in
   `csa/recon/diffusion.py`; wire it to the csgm-mri-langevin sampler using
   the checkpoint above.
6. **ESPIRiT** via `sigpy` for real-data coil maps (`csa/physics/espirit.py`
   is written but untested against real fastMRI files, since sigpy is not
   installed).
7. **fastMRI+ annotations**: done — cloned 2026-09-08 to
   `/Volumes/T9/fastMRI_brain/fastmri_plus/Annotations/brain.csv`.
8. **The pre-registration and thresholds** in `configs/thresholds.yaml`, still
   null (see `prereg/README.md`) — required before any pair is scored, not
   before code can run.
