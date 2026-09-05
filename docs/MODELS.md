# Models under audit

| Family | Implementation | Weights | Output | Status |
|---|---|---|---|---|
| E2E-VarNet | `fastmri.models.VarNet` | public fastMRI brain leaderboard checkpoint | magnitude (RSS) by default; the audit needs the complex cascade output, so the adapter must expose it | adapter written, untested |
| Complex U-Net | `fastmri.models.Unet(in_chans=2, out_chans=2)` | trained here on the standard fastMRI brain split | complex | training script to write |
| Diffusion | csgm-mri-langevin (Jalal et al 2021) | released brain multi-coil weights | complex | adapter stub |

Rules: models are audited as released; no model is retrained on the lesion
bank; no lesion-bearing slice enters any training or calibration step; both
arms of a pair use the same seed; seeds are then varied across replicates.

Chung and Ye's score-MRI weights were trained on knee data and are therefore
not the first choice for a brain audit.
