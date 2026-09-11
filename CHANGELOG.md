# Changelog

All notable changes to this repository. Dates are ISO. Every deviation from the
frozen pre-registration must be recorded here with its reason.

## Unreleased

### Amendments A2 and A3 (2026-09-12, authors' decision)

- **A2, H1's comparator.** The endpoint stays the silent erasure rate; the
  test becomes a two-sided paired contrast against the no-prior baseline
  (`cg_sense_tuned`) on the same lesions, acquisitions and seeds, reported
  overall and at matched `mu_lambda`. `pi0` becomes a descriptive reference.
  Reason: the classical audit showed a threshold test against `pi0` is
  satisfied by the acquisition alone. Sizing targets a half-width of 0.05 on
  the contrast, which is not more demanding than the original interval on a
  single rate because the paired contrast of correlated rates has smaller
  variance.
- **A3, scope.** All claims are stated as being about focal structure of
  specified size and contrast under a specified acquisition; the clinical
  findings that motivate the work are named as motivation only. The bank's
  conspicuous contrast levels are about four times the measured core contrast
  of the fastMRI+ lacunar-infarct boxes, which are predominantly chronic, and
  that extrapolation is now stated in the limitations.
- Manuscript updated accordingly (Introduction, causal graph label, Case A
  assumption, H1, statistical plan, limitations, conclusion note); 36 pages.
  Drafts in `prereg/addendum_v1_DRAFT.md`; `configs/confirmatory_A1b.yaml`
  carries the H1 specification and the baseline arm.


### Fixed (2026-09-11) — measurement correction, affects reported silent erasure

Global PSNR and SSIM scored each arm of a pair against its own reference, so a
bright lesion raised the reference maximum and shifted the score by
20 log10(max_cf / max_f) with no change in the error (worst case seen: 0.871 dB
recorded against 0.0006 dB of real change). Both arms are now scored with
`maxval` fixed to the lesion-absent reference; `tests/test_common_scale.py`
pins it. Erasure rates are unaffected (observer only); silent erasure was
understated. The affected table is quarantined under
`outputs/annotated_val_A1/superseded_buggy_psnr/`.

### Added

- `cg_sense_tuned`: CG-SENSE at an image-quality-oriented lambda, the fair
  no-prior baseline (the audit's noise-set lambda amplifies noise by ~10 dB
  and inflates its erasure rate).


### Amendment A1 to the lesion bank (2026-09-10) — DEVIATION from pre-registration v1

Reason: the pre-registered bank (1 to 27 mm^3 at fastMRI+-derived contrasts)
is 97 % undetectable in the fully sampled reference on 5 to 7.5 mm slices, so
the erasure endpoint would be empty (docs/VALIDATION.md, 2026-09-09; the
classical dry run confirmed 0 eligible pairs of 714). Change: volumes 10
(sub-threshold control), 27, 64, 100, 200 mm^3; FLAIR contrasts 0.24 (fastMRI+
focal median), 1.0, 1.5; T1 contrasts -0.21, -0.6, -1.0; insertion sites
restricted to the ESPIRiT support. Unchanged: every threshold, the erasure
definitions, the observer, the hypotheses and the analysis plan. Draft
addendum for the authors in `prereg/addendum_v1_DRAFT.md`; not deposited.
Verified on the reference 2026-09-10: 25 % of planned insertions eligible
(FLAIR 1.5 majority-eligible from 27 mm^3; T1 capped near 44 %); option A1a
(leaner confirmatory grid) proposed in the draft.
Decided by the authors ("do as you suggest", 2026-09-10) on the basis of the
reference-detectability grid; 44 lacunar-infarct boxes measured from 13
fastMRI+ files were found to be chronic, modest-contrast lesions and so could
not supply acute-like levels on their own.


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

### Overnight CPU stage completed (2026-09-09 06:15)

- Characterization of all 984 planned insertions, Experiment 2 and 3
  demonstrations on the annotated set, validation figures
  (`scripts/09_figures.py`). Headline: only 3 % of the planned bank is
  detectable in the fully sampled reference under the pre-registered
  `z_det`, because 1 to 27 mm^3 lesions are diluted by 5 to 7.5 mm slices.
  Recorded in `docs/VALIDATION.md` for the sizing addendum; no threshold or
  definition changed. `scripts/explore_reference_detectability.py` added to
  quantify the options (larger volumes, higher contrasts, thin-slice
  equivalent).

### Experiment 1 classical dry run completed (2026-09-09 09:20)

- 714 pairs, zero-filled and CG-SENSE at R = 8, on 15 annotated volumes.
  PSNR changes by at most 0.023 dB from the lesion and stays inside the
  lesion-free interval in 94 % of pairs; classical reconstruction retains
  0 to 69 % of reference detectability; no pair passes the reference gate.
  Validation entry, not a result. Dry-run figure added to `09_figures.py`.
- `02_bank.py` restricts insertion sites to the ESPIRiT support when a cached
  map exists (18 sites fell outside it in the dry run).

### Changed (2026-09-10, hand-off)

- Global PSNR/SSIM in `04_exp1.py` on the fastMRI 320 x 320 centre crop by
  default (`--metric-crop`), matching how the field reports them.
- `scripts/train_unet_complex.py` implemented and CPU smoke-tested;
  `csa/recon/unet.py` reads the architecture from the checkpoint.

### Changed

- `scripts/04_exp1.py` implemented: Algorithm 1 per pair for any reconstructor,
  cross-fitted observers on reconstructions, erasure flags from the
  pre-registered thresholds, resumable pass 1, `--observer-only` pass 2. Run
  labels keep the classical CPU dry run apart from the confirmatory GPU run.
  Smoke-tested on real AXT2 slices; under load one (file, slice, model, R) unit
  takes about 26 minutes on CPU, so the overnight dry run is R = 8 only, two
  slices per file, zero-filled and CG-SENSE.
- `scripts/02_bank.py`: lesion-core contrast statistic over focal fastMRI+
  labels replaces the box mean, which loose boxes diluted to zero. FLAIR
  levels [0.12, 0.24, 0.36] from 14 focal boxes (median 0.24, quartiles
  collapsed so widened to x0.5 / x1.5); T1 has 3 focal boxes and uses the
  signed defaults, recorded as `default_insufficient_focal_boxes`.
- Overnight chain on the 17 annotated val volumes: integrity -> bank ->
  characterization (82 slice units, 4 workers) -> Experiment 2 and 3
  demonstrations -> Experiment 1 classical dry run.
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
