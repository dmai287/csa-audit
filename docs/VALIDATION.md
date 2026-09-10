# Validation log

Checks of the audit machinery that are **not experimental results**. Nothing
here applies the pre-registered thresholds, scores an erasure endpoint, or
touches the output of a learned reconstructor. Entries are dated; the code
state is the commit named in each entry.

## 2026-09-08: statistics machinery, by simulation (`scripts/validate_stats.py`)

Design under test: `pi0 = 0.05`, intra-subject correlation `rho = 0.10`,
`m = 12` pairs per subject; the pairs-per-cell rule returns 958 pairs from 80
subjects. 100 replications, 300 bootstrap resamples each.

| Check | Result |
| --- | --- |
| Cluster-bootstrap 95% interval, coverage of the true rate | 0.95 |
| Naive binomial interval ignoring clustering, coverage | 0.81 |
| Empirical half-width from the sample-size rule | 0.0198 (target <= 0.02) |
| Equivalence rule (TOST, margin 0.2), true SMD 0.5 | declared in 0 of 50 runs |
| Equivalence rule (TOST, margin 0.2), true SMD 0.0, n = 400 | declared in 17 of 50 runs |

Reading. The interval and the sizing rule behave as the analysis plan claims,
and the naive interval's under-coverage is the reason the plan resamples
subjects. The last row is a power statement, not a defect: with 400 pairs the
90% interval for a standardised mean difference is about +/-0.17 before
cluster inflation, so equivalence within +/-0.2 is declared only when the
estimate falls near zero. Two consequences for the sizing addendum:

1. H2 must be tested on the pooled confirmatory sample, not per cell.
2. H2's power depends on the number of *erased* pairs, which is the erasure
   rate times the sample, not on the sample alone. At `pi0 = 0.05` and 958
   pairs per cell pooled over nine cells, the erased group is a few hundred
   pairs and the standardised-difference standard error is roughly 0.05,
   which gives adequate power at margin 0.2. If the realised erasure rate is
   much lower than `pi0`, H2 becomes under-powered even when H1 is not, and
   the addendum should say so before the confirmatory run.

## 2026-09-08: audit physics on a real acquisition (commit 85de9e8)

`file_brain_AXT2_202_2020533.h5`, slice 8, 12 coils, 640x320, equispaced
R = 8, classical CG-SENSE only.

| Check | Result |
| --- | --- |
| ESPIRiT maps (sigpy), peak `sum|S|^2` / support / time | 1.000 / 28% of FOV / 38 s |
| Adjoint identity of the real operator | relative error 4.4e-15 |
| Two insertion paths (Algorithm 1) | max difference 0 |
| Noise cancellation between arms | max residual 5.5e-19 |
| Measurement gain, 10 mm^3 lesion | `kappa^2` = 0.235 (single-coil 1/R = 0.125) |
| Measurement fraction | `mu_lambda` = 0.426 at `lambda` = 1.4e-3 |
| CG-SENSE lesion transfer | `t_R` = 0.998, `t_N` = 0.032 |

Reading. The measurement-gain result supports the manuscript's Section 3.4
claim that coil encoding raises a compact lesion's measured fraction above
`1/R` without restoring it to one. `t_R` near 1 with `t_N` near 0 is the
expected signature of a linear data-consistent reconstruction, which supplies
almost nothing in the unmeasured subspace.

Timing on CPU: about 4 minutes per pair for the classical reconstructor and
the null-space statistics alone. `mu_lambda` and `P_lambda(ell)` do not depend
on the model and are computed once per lesion by
`scripts/02b_characterize_bank.py`.

## 2026-09-08: fastMRI+ brain annotation coverage

From `Annotations/brain.csv` (8,213 rows, 997 files): sequences AXFLAIR
(4,210 rows), AXT1 (2,322), AXT1POST (1,681); no AXT2; 113 rows labelled
"Lacunar infarct"; no hemorrhage label. Consequence: T2 hosts carry the
synthetic (Case A) bank only, with default contrast levels recorded as
`contrast_source = default`, and have no Case B counterpart.

## 2026-09-08: measurement-side acquisition shift (`scripts/exp3_measurement_side.py`)

VALIDATION, NOT EXPERIMENT 3'S RESULT. Four lesions (5, 10, 10, 27 mm^3) on
two AXT2 slices, R = 8, twelve interventions on the acquisition at fixed
anatomy and lesion; means over lesions. `kappa^2` is the measurement gain,
`mu_lambda` the measured fraction under the noise-set filter.

| Condition | `kappa^2` | `mu_lambda` | sampled fraction |
| --- | --- | --- | --- |
| `equispaced_offset_0` | 0.236 | 0.486 | 0.128 |
| `equispaced_offset_1` | 0.236 | 0.486 | 0.125 |
| `equispaced_offset_2` | 0.236 | 0.482 | 0.125 |
| `equispaced_offset_3` | 0.236 | 0.475 | 0.125 |
| `random_seed_0` | 0.240 | 0.422 | 0.128 |
| `random_seed_1` | 0.215 | 0.407 | 0.134 |
| `random_seed_2` | 0.239 | 0.467 | 0.125 |
| `noise_x2` | 0.236 | 0.416 | 0.128 |
| `noise_x4` | 0.236 | 0.351 | 0.128 |
| `coils_first_half` | 0.236 | 0.396 | 0.128 |
| `coils_second_half` | 0.236 | 0.493 | 0.128 |
| `coils_alternate` | 0.236 | 0.448 | 0.128 |

Reading.
- Shifting the equispaced pattern by one to three lines leaves both statistics
  essentially unchanged (mu 0.475 to 0.486).
- Random masks at the same sampled fraction measure the lesion less and less
  consistently (mu 0.41 to 0.47 across seeds).
- Raising the noise level lowers `mu_lambda` (0.486 to 0.416 to 0.351 at x2 and
  x4) while `kappa^2` is untouched: the gain is a property of the operator, the
  measured fraction is a property of the operator and the noise, which is
  what the lambda rule is for.
- Coil subsets change `mu_lambda` (0.40 to 0.49) but not `kappa^2`, because the
  subset maps are renormalised to `sum|S|^2 = 1`: with normalised maps the gain
  depends on the mask and the lesion spectrum only, and coil geometry enters
  through the conditioning of the operator, i.e. through `mu_lambda`. This is
  worth a sentence in the manuscript when the two statistics are introduced.

Cost: about 35 s per condition per lesion on CPU (one CG solve each).

## 2026-09-08: residual insensitivity on real acquisitions (`scripts/exp2_physics_demo.py`)

VALIDATION, NOT EXPERIMENT 2'S RESULT. Eight lesions (1 to 27 mm^3) on two
AXT2 slices, R = 8 equispaced, classical CG-SENSE at the noise-set lambda;
reference known, so the reconstruction error `e` can be split into its
filtered range and null-space parts.

| Quantity (mean over 8 lesions) | Value |
| --- | --- |
| Share of error energy in the null space | 0.55 |
| Relative residual of the reconstruction | 0.0822 |
| Residual after removing the null-space error component | 0.0828 |
| Residual after removing the range-space error component | 0.0868 |
| Pure null-space edit: change in residual | 2.4e-08 |
| Pure null-space edit: Proposition 1 bound `||A delta||/||Y||` | 1.2e-05 |
| Pure null-space edit: `||A delta||/||delta||` | 0.010 |
| CG-SENSE transfer `t_R`, `t_N` | 0.997, 0.042 |
| Proposition 1 bound held | 8 of 8 |

Reading. More than half of the classical reconstruction's error lives in the
unmeasured subspace, and the k-space residual does not see it: deleting that
component changes the residual by under one percent, while deleting the
measured component moves the residual to the noise floor (the reconstruction
had fitted part of the noise). A lesion-shaped edit confined to the numerical
null space passes through `A` at about one percent of its norm and moves the
residual by about 1e-8, three orders of magnitude inside the bound. This is
Section 3.4 of the manuscript exhibited on measured data with a linear
reconstructor; Experiment 2 proper repeats the decomposition on the learned
models' erased pairs.

Cost: about 3 minutes per lesion on CPU (two CG-SENSE reconstructions and
three filtered-projector solves).

## 2026-09-09: overnight CPU stage on 17 annotated val volumes (run `annotated_val0`)

Integrity (17 volumes, none excluded) -> bank -> characterization of 984
planned insertions (328 site-volume lesions x 3 accelerations; 82 slice units,
4 workers, 4 h 05 min) -> Experiment 2 and 3 demonstrations. Outputs in
`outputs/annotated_val0/`, figures in `outputs/annotated_val0/figures/`.

### Model-independent measurement of the planned lesions

Mean over lesions; the single-coil prediction is `1/R` = 0.250, 0.167, 0.125.

| Volume | `kappa^2` R=4 / 6 / 8 | `mu_lambda` R=4 / 6 / 8 |
| --- | --- | --- |
| 1 mm^3 | 0.32 / 0.23 / 0.17 | 0.95 / 0.74 / 0.59 |
| 5 mm^3 | 0.38 / 0.27 / 0.20 | 0.96 / 0.75 / 0.62 |
| 10 mm^3 | 0.44 / 0.32 / 0.23 | 0.97 / 0.78 / 0.63 |
| 27 mm^3 | 0.49 / 0.37 / 0.27 | 0.96 / 0.79 / 0.67 |

Both statistics are monotone in volume and in acceleration, identical between
FLAIR and T1 hosts to two decimals, and the gain sits 1.3x (1 mm^3) to 2.2x
(27 mm^3, R=8) above the single-coil `1/R` without approaching one: the
Section 3.4 claim, now measured on 328 real-acquisition lesions. CG note:
median 80 iterations (the cap) with worst relative residual 1.2e-3; raise
`--cg-maxiter` on the GPU run.

### Reference detectability of the planned bank (the finding that matters)

Cross-fitted CHO on the fully sampled reference, `z_det` = 3.0 as
pre-registered. Only **3 %** of planned lesions reach `z_det` (1 mm^3: 0 %,
5 mm^3: 4 %, 10 mm^3: 2 %, 27 mm^3: 5 %). Median `z_ref` is 0.1 to 0.6 for
FLAIR at contrasts 0.12 to 0.36 and 0.5 to 2.1 for T1 at contrasts -0.3 to
-1.0; reference `d'` is below 0.7 for every FLAIR cell and below 2.2 for T1.

Why: fastMRI brain slices are 5 to 7.5 mm thick, so a sphere of 10 mm^3
(radius 1.3 mm) fills 0.24 of the slice and a 27 mm^3 sphere 0.33; the
effective in-plane contrast of a 0.24 FLAIR lesion is about 0.06. The
manuscript's own motivating example is a 100 mm^3 lacunar infarct; the bank's
1 to 27 mm^3 range sits below what this data can show even without
undersampling. Under the pre-registered erasure definition (z_ref >= z_det
required) the current bank would yield almost no eligible pairs.

This is a pilot-stage finding to be resolved by the authors in the
pre-registration's sizing addendum before any confirmatory run, not by the
code. Options being quantified by `scripts/explore_reference_detectability.py`:
extend the volume range upward (64 to 500 mm^3); higher contrast levels
(lacunar-infarct boxes exist in fastMRI+ but not in these 17 volumes); or a
site-matched observer, which would be a definitional change and must be
declared as such. The thresholds themselves were not touched.

### Experiment 2 demonstration on the annotated set (24 lesions, R=8, CG-SENSE)

Bound held 24 of 24. Null-space share of the error 0.30 (0.21 to 0.52);
deleting it moves the residual by 0.6 %, deleting the measured part moves it
to the noise floor; a null-space lesion edit passes through `A` at 0.6 % of its
norm and does not change the residual at the reported precision. `t_R` 0.999,
`t_N` 0.08 (0.03 to 0.27).

### Experiment 3 measurement side on the annotated set (12 lesions, R=8)

`kappa^2` 0.225 under every condition except random masks; `mu_lambda` 0.61
to 0.63 for equispaced offsets, noise x2 and x4, and coil subsets, but 0.49
to 0.58 for random masks. On these volumes the noise interventions barely
moved `mu_lambda` (they did on the AXT2 smoke set), which says the noise-set
`lambda` is small relative to the measured singular values here.

## 2026-09-09: reference detectability over an extended grid (`scripts/explore_reference_detectability.py`)

Cross-fitted CHO on the fully sampled reference, 82 slices, 328 sites, 4,428 patches; cells show
median `z_ref` and, in brackets, the fraction of lesions at or above the pre-registered `z_det` = 3.0.
'thick slice' applies the partial-volume fill of the actual 5-7.5 mm slices; 'thin-slice equivalent' does not.

**AXFLAIR, thick slice (as acquired)**

| contrast | 10 mm^3 | 27 mm^3 | 64 mm^3 | 100 mm^3 | 200 mm^3 | 500 mm^3 |
|---|---|---|---|---|---|---|
| +0.24 | 0.3 (1%) | 0.4 (2%) | 0.5 (3%) | 0.5 (3%) | 0.6 (4%) | 1.0 (4%) |
| +0.50 | 0.7 (4%) | 0.9 (5%) | 1.1 (7%) | 1.3 (8%) | 1.5 (15%) | 2.4 (28%) |
| +1.00 | 1.7 (14%) | 2.0 (22%) | 2.5 (37%) | 2.9 (47%) | 3.3 (59%) | 4.9 (87%) |
| +1.50 | 2.6 (36%) | 3.2 (58%) | 4.1 (75%) | 4.4 (77%) | 5.2 (84%) | 7.5 (94%) |

**AXFLAIR, thin-slice equivalent**

| contrast | 10 mm^3 | 27 mm^3 | 64 mm^3 | 100 mm^3 | 200 mm^3 | 500 mm^3 |
|---|---|---|---|---|---|---|
| +0.24 | 1.2 (3%) | 1.1 (4%) | 1.1 (4%) | 1.0 (5%) | 0.9 (6%) | 1.1 (4%) |
| +0.50 | 2.4 (30%) | 2.3 (23%) | 2.2 (24%) | 2.2 (26%) | 2.1 (27%) | 2.5 (32%) |
| +1.00 | 5.1 (92%) | 4.8 (87%) | 4.7 (86%) | 4.7 (83%) | 4.5 (83%) | 5.1 (89%) |
| +1.50 | 7.9 (98%) | 7.5 (96%) | 7.3 (96%) | 7.1 (94%) | 7.1 (91%) | 7.9 (96%) |

**AXT1, thick slice (as acquired)**

| contrast | 10 mm^3 | 27 mm^3 | 64 mm^3 | 100 mm^3 | 200 mm^3 | 500 mm^3 |
|---|---|---|---|---|---|---|
| -0.30 | 0.8 (0%) | 0.9 (0%) | 0.9 (1%) | 1.0 (0%) | 1.1 (1%) | 1.7 (2%) |
| -0.60 | 1.5 (4%) | 1.7 (7%) | 1.7 (6%) | 1.7 (8%) | 2.1 (16%) | 3.2 (55%) |
| -1.00 | 2.3 (25%) | 2.8 (41%) | 2.8 (40%) | 2.8 (37%) | 3.4 (56%) | 5.2 (95%) |

**AXT1, thin-slice equivalent**

| contrast | 10 mm^3 | 27 mm^3 | 64 mm^3 | 100 mm^3 | 200 mm^3 | 500 mm^3 |
|---|---|---|---|---|---|---|
| -0.30 | 1.8 (1%) | 1.7 (0%) | 1.5 (0%) | 1.4 (0%) | 1.4 (1%) | 1.7 (2%) |
| -0.60 | 3.3 (68%) | 3.1 (60%) | 2.7 (36%) | 2.4 (26%) | 2.5 (26%) | 3.2 (63%) |
| -1.00 | 5.3 (97%) | 5.1 (96%) | 4.6 (89%) | 4.1 (84%) | 4.1 (83%) | 5.4 (97%) |

Reading.
- Contrast, not volume, governs reference detectability on this data. At the fastMRI+-derived FLAIR
  contrast of 0.24 no volume up to 500 mm^3 is detectable on thick slices (<= 6 %); at contrast 1.0
  a majority of lesions of 200 mm^3 and above are, and at 1.5 a majority of 64 mm^3 and above.
- Partial volume is the mechanism: the thin-slice equivalent at contrast 1.0 puts 83 to 92 % of FLAIR
  lesions of every volume above `z_det`, and 83 to 97 % of T1 lesions at contrast -1.0.
- The manuscript's motivating 100 mm^3 lacunar infarct on a 7.5 mm slice (fill 0.51) reaches
  `z_ref` 2.9 at contrast 1.0 and 4.5 at 1.5.
- Without dilution `z_ref` is nearly flat in volume at fixed contrast: with anatomy-dominated
  background covariance, a larger lesion also sees more background variability at its own scale.

Implication for the sizing addendum (the authors' decision): the 0.24 FLAIR level came from
'nonspecific white matter lesion' boxes, not from the acute-stroke lesions the paper is about;
fastMRI+ holds 113 'Lacunar infarct' boxes, none in these 17 volumes. Deriving the contrast levels
from those boxes, and moving the bank's volumes to 27-200 mm^3 with 10 mm^3 retained as a
sub-threshold control, would give a bank whose majority is detectable in the reference. Both are
bank-parameter choices the pre-registration leaves to the addendum; `z_det` itself is unchanged.

## 2026-09-09: Experiment 1 machinery, classical reconstructors (run `dryrun_classical`)

DRY RUN, NOT A CONFIRMATORY RESULT. 714 pairs (119 lesions x 3 contrasts x
2 reconstructors) on 15 annotated volumes, two central slices each, R = 8
equispaced; both arms reconstructed, observers cross-fitted on the
reconstructions, pre-registered erasure logic applied. 60 units, 2 h 14 min
on 4 CPU workers; 18 sites skipped for lying outside the coil-map support.

| | zero-filled | CG-SENSE |
| --- | --- | --- |
| `t_R` (measured transfer) | 0.769 | 0.999 |
| `t_N` (null-space transfer) | 0.000 | 0.096 |
| relative residual | 0.155 | 0.125 |
| PSNR, lesion-absent arm | 27.4 dB | 17.5 dB |
| lesion-present PSNR inside the lesion-free 95 % interval | 94 % of pairs | 94 % of pairs |
| median / max PSNR change from the lesion | 0.000 / 0.020 dB | 0.000 / 0.023 dB |
| detectability retained, median `z_recon / z_ref` where `z_ref` > 1, FLAIR / T1 | -0.07 / 0.69 | 0.27 / 0.52 |
| ensemble `d'`, T1 27 mm^3 at contrast -1.0: reference 2.16 | 1.02 | 0.28 |
| pairs passing the reference gate `z_ref >= z_det` | 0 of 357 | 0 of 357 |

Reading.
- The global score is blind to the lesion on real data: across 714 pairs the
  lesion changes PSNR by at most 0.023 dB, the order predicted by the
  companion review's identity for lesions of this size, and it stays inside
  the lesion-free interval in 94 % of pairs for both reconstructors.
- Both classical reconstructions lose most of the reference detectability at
  R = 8, in the direction the theory predicts: zero-filling supplies nothing
  in the null space (`t_N` = 0) and its measured transfer is 0.77 because the
  adjoint is not data-consistent; CG-SENSE writes the measured part back
  (`t_R` = 0.999) but at the noise-set lambda amplifies noise, which is why its
  PSNR is 10 dB lower and its ROI contrast collapses.
- No pair passes the reference gate, so the pre-registered erasure endpoint is
  empty for the current bank: the amendment case, seen end to end.

Figures: `outputs/annotated_val0/figures/dryrun_exp1_classical.png`.

## 2026-09-10: amended bank (A1) verified on the reference

Patch-only check on 115 slices of 22 hosts (2,542 patches), cross-fitted CHO,
pre-registered `z_det` = 3.0. Cells: fraction of planned lesions at or above
`z_det` (median `z_ref` in brackets), thick slices as acquired.

| sequence, contrast | 10 mm^3 | 27 mm^3 | 64 mm^3 | 100 mm^3 | 200 mm^3 |
|---|---|---|---|---|---|
| AXFLAIR +0.24 | 0 % (z 0.3) | 3 % (z 0.4) | 5 % (z 0.5) | 3 % (z 0.5) | 3 % (z 0.7) |
| AXFLAIR +1.00 | 15 % (z 1.6) | 26 % (z 2.1) | 35 % (z 2.6) | 46 % (z 2.9) | 63 % (z 3.4) |
| AXFLAIR +1.50 | 38 % (z 2.6) | 57 % (z 3.2) | 76 % (z 4.0) | 79 % (z 4.5) | 85 % (z 5.1) |
| AXT1 -0.21 | 2 % (z 0.6) | 1 % (z 0.8) | 1 % (z 0.8) | 0 % (z 0.7) | 0 % (z 0.8) |
| AXT1 -0.60 | 4 % (z 1.4) | 5 % (z 1.7) | 5 % (z 1.8) | 3 % (z 1.6) | 7 % (z 1.7) |
| AXT1 -1.00 | 27 % (z 2.3) | 44 % (z 2.8) | 41 % (z 2.8) | 30 % (z 2.6) | 44 % (z 2.8) |

Expected eligible fraction of planned insertions: 25 % overall
(FLAIR 36 %, T1 14 %), about 429 of 1725,
up from 3 % for the pre-registered bank. Reading: FLAIR at 1.5 is majority-
eligible from 27 mm^3 up and FLAIR at 1.0 only at 200 mm^3; T1 is capped
near 44 % at full hypointensity, a limit of thick slices rather than a
parameter choice; the subtle strata (FLAIR 0.24, T1 -0.21) are sub-threshold
as intended. Consequence for sizing: the erasure endpoint's denominator is
the eligible pairs, so the pairs-per-cell rule applies to eligible pairs and
the confirmatory bank should oversample the eligible cells. A leaner
confirmatory grid is proposed in the addendum draft as option A1a.

## 2026-09-10: model adapters checked on CPU before hand-off

- **Complex U-Net training script**: one epoch on two real slices with
  validation and checkpointing (CPU, 140 s); the adapter rebuilt the
  architecture from the checkpoint and returned a finite complex image on a
  real 12-coil slice. The real model is trained on the GPU machine.
- **Global metrics convention**: fastMRI reports PSNR and SSIM on the central
  320 x 320 crop, not the full 640 x 320 grid. On the full grid a
  reconstruction that fills the empty readout margins differently from the
  support-masked reference is penalised for nothing diagnostic, which the
  VarNet check exposed. `scripts/04_exp1.py` now scores global PSNR/SSIM on
  that crop by default (`--metric-crop 320`; 0 restores the full grid). The
  classical dry run used the full grid; its PSNR-change result (max 0.023 dB)
  is unaffected in kind, since the lesion sits inside the crop.
