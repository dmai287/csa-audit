# Changelog

All notable changes to this repository. Dates are ISO. Every deviation from the
frozen pre-registration must be recorded here with its reason.

## 0.1.0 — 2026-09-05

- Initial scaffold. NumPy core: multi-coil SENSE operator with adjoint test;
  fastMRI-convention equispaced and random masks; additive lesion insertion
  (Algorithm 1) with two-path equality test; conjugate-gradient filtered
  range/null-space projectors with closed-form and transfer-identity tests;
  channelised Hotelling observer with Laguerre-Gauss channels and an
  ideal-observer bound test; PSNR/SSIM in fastMRI convention, ROI CNR and ROI
  SSIM; pre-specified erasure and silent-erasure logic; cluster bootstrap,
  two-one-sided-tests equivalence, AUC intervals and the pairs-per-cell rule.
- Reconstructor adapters: zero-filled and CG-SENSE (NumPy); E2E-VarNet,
  complex U-Net and Langevin diffusion adapters behind lazy torch imports,
  untested against real weights.
- Scripts: synthetic end-to-end smoke run, fastMRI integrity check, pilot
  sample-size calculator, experiment and statistics skeletons.
- No experiment has been run. No result exists.
