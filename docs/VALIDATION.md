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
