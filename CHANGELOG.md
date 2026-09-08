# Changelog

All notable changes to this repository. Dates are ISO. Every deviation from the
frozen pre-registration must be recorded here with its reason.

## Unreleased

### First real-data validation (2026-09-08)

The audit physics ran end to end on a genuine fastMRI brain multi-coil slice
(`file_brain_AXT2_202_2020533.h5`, slice 8, 12 coils, 640x320, R=8 equispaced),
with no reconstruction model involved and no experimental endpoint computed.
Results, all consistent with the theory in the manuscript:

| Check | Result |
| --- | --- |
| ESPIRiT maps (sigpy) | `sum|S|^2` peak 1.000, support 28% of FOV, 38 s/slice |
| Adjoint identity on the real operator | relative error 4.4e-15 |
| Two lesion-insertion paths (Algorithm 1) | agree exactly, max difference 0 |
| Noise cancellation between arms | max residual 5.5e-19 |
| Lesion measurement gain, 10 mm^3 at R=8 | `kappa^2` = 0.235 against 1/R = 0.125 |
| Lesion measurement fraction | `mu_lambda` = 0.426 |
| CG-SENSE lesion transfer | `t_R` = 0.998, `t_N` = 0.032 |

The measurement-gain result is direct support for the manuscript's claim that
coil encoding raises a compact lesion's measured fraction above the
single-coil `1/R` prediction without restoring it to one. `t_R` near 1 with
`t_N` near zero is the expected signature of a linear data-consistent
reconstruction: it writes back what was measured and supplies almost nothing
in the unmeasured subspace, which is why `t_N` is the statistic that
distinguishes learned priors.

Timing, measured on CPU (this machine has no CUDA device): about 4 minutes per
pair for the classical reconstructor alone. Two consequences for the
confirmatory runs, both recorded here rather than discovered later:
a GPU is required at experiment scale, and `P_lambda(ell)` / `mu_lambda`
depend only on the lesion and the operator, not on the model, so they should
be computed once per (lesion, acquisition) and reused across models.

### Validation on real data, continued (2026-09-08, evening)

- Experiment 2's operator-level content and Experiment 3's measurement-side
  content demonstrated on real AXT2 slices with CG-SENSE; statistics machinery
  checked by simulation. All in `docs/VALIDATION.md`. None of it applies the
  pre-registered thresholds or scores an endpoint.
- Selective extraction of 30 annotated FLAIR/T1 volumes from val batch 0 and a
  detached chain that runs the full CPU stage (`scripts/run_cpu_stage.sh`) on
  them when extraction completes.

### Changed

- `csa/physics/espirit.py` now rejects input that is not `(coils, H, W)`. A
  4-D array was silently treated by sigpy as a single-coil 3-D volume and
  returned all-zero maps, which propagated NaNs through every downstream
  statistic. Caught during the validation above.
- `scripts/02_bank.py` implemented against the real fastMRI+ schema
  (`file, slice, study_level, x, y, width, height, label`): measures per-box
  appearance statistics from real data, derives per-sequence contrast levels
  from their interquartile range, and samples insertion sites from
  annotation-free slices. fastMRI+ brain covers AXFLAIR, AXT1 and AXT1POST
  only, with no AXT2 and no hemorrhage label; 113 boxes are labelled
  "Lacunar infarct". T2 hosts therefore support the synthetic (Case A) bank
  but have no Case B real-lesion counterpart.
- `scripts/train_unet_complex.py` added as a scoped stub (needs torch, fastmri
  and CUDA).

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
