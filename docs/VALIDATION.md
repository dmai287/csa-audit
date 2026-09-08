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
