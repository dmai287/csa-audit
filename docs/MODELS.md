# Models under audit

| Family | Implementation | Weights | Output | Status |
|---|---|---|---|---|
| E2E-VarNet | `fastmri.models.VarNet` | downloaded 2026-09-08 to `/Volumes/T9/fastMRI_brain/checkpoints/varnet_brain_leaderboard_state_dict.pt` (114 MB, from `dl.fbaipublicfiles.com`, the official fastMRI leaderboard checkpoint; the fastmri.org leaderboard site itself is down post Meta-to-NYU handover, but this direct file host still serves the weights) | magnitude (RSS) by default; the audit needs the complex cascade output, so the adapter must expose it | weights in hand; adapter still needs the complex-output change before use |
| Complex U-Net | `fastmri.models.Unet(in_chans=2, out_chans=2)` | none published; train on `train_batch_0` (already downloaded) | complex | training script to write (`scripts/train_unet_complex.py`) |
| Diffusion | csgm-mri-langevin (Jalal et al 2021) | downloaded 2026-09-08 to `/Volumes/T9/fastMRI_brain/checkpoints/ncsnv2-mri-mvue/` (1.3 GB, from the repo's published Google Drive link); one unconditional NCSNv2 score model, used across contrasts via sampling config, not per-contrast checkpoints | complex | weights in hand; adapter is still a stub (`csa/recon/diffusion.py`) and the repo's sampler needs BART installed for sensitivity maps |

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
