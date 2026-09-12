# Sizing addendum v1 — DRAFT for author review (not deposited)

Prepared 2026-09-10 from the CPU-only stage of the audit. This document is
a draft for the authors to edit and deposit on OSF as the addendum that
pre-registration v1 anticipates ("the realised pairs-per-cell is posted as a
dated OSF addendum before any confirmatory reconstruction"). Nothing in it
has been deposited, and no confirmatory reconstruction has been run.

## What was found before any reconstruction

The planned bank was characterised on the fully sampled reference alone
(`scripts/02b_characterize_bank.py`; `docs/VALIDATION.md`, entries of
2026-09-09). Two facts follow from the data, not from any model:

1. **The pre-registered bank is largely undetectable in the reference.** With
   volumes 1, 5, 10 and 27 mm^3 at contrast levels derived from fastMRI+
   boxes (FLAIR 0.12 to 0.36; T1 -0.3 to -1.0), only 3 % of planned lesions
   reach the pre-registered `z_det` = 3.0 in the fully sampled reference.
   Under the pre-registered erasure definition, which requires
   `z_ref >= z_det`, the Experiment 1 endpoint would be empty; the classical
   dry run (714 pairs) confirmed 0 eligible pairs.
2. **The cause is partial volume on 5 to 7.5 mm slices**, compounded by the
   modest contrast of the lesions annotated in fastMRI+. A 10 mm^3 sphere
   fills 0.24 of its slice; the thin-slice equivalent of the same bank puts
   83 to 97 % of lesions above `z_det`. The fastMRI+ "Lacunar infarct" boxes
   (44 measured in 13 files) are chronic lacunes with modest, mixed-sign core
   contrast on FLAIR (median -0.22, range -0.63 to +0.97), not the bright
   acute infarcts the audit targets.

## Amendment A1 to the bank (a bank-parameter change, not a threshold change)

| Parameter | Pre-registered | Amended |
| --- | --- | --- |
| Volumes (mm^3) | 1, 5, 10, 27 | 10 (sub-threshold control), 27, 64, 100, 200 |
| FLAIR contrast levels | fastMRI+ focal IQR (0.12, 0.24, 0.36) | 0.24 (fastMRI+ focal median, "subtle"), 1.0, 1.5 ("acute-like") |
| T1 contrast levels | defaults (-0.3, -0.6, -1.0) | -0.21 (fastMRI+ focal median), -0.6, -1.0 |
| Insertion sites | RSS-thresholded tissue | additionally restricted to the ESPIRiT map support |

Justification for the acute-like levels: the reference-detectability grid
(`docs/VALIDATION.md`, 2026-09-09) shows that on these slices a majority of
lesions of 64 mm^3 and above clear `z_det` at contrast 1.0 to 1.5, whereas no
volume does at 0.24. Acute ischaemic lesions are markedly brighter on FLAIR
and DWI than the chronic lesions annotated in fastMRI+, so the amended axis
spans subtle-chronic to acute-like appearance rather than the chronic range
alone. The subtle level and the 10 mm^3 volume are retained as strata that
are expected to fail the reference gate, so that the audit also records what
happens to lesions near and below detectability.

Unchanged: `z_det` = 3.0, `z_miss` = 2.0, `z_fs` = 4.0, `pi0` = 0.05, the
erasure and silent-erasure definitions, the observer (cross-fitted CHO with
Laguerre-Gauss channels), the `lambda` rule, the hypotheses, the equivalence
margins, the multiplicity rule and the pilot design.

## Verification of the amended bank on the reference (to be filled)

`outputs/annotated_val_A1/reference_detectability_check.csv`, 115 slices of
22 hosts, cross-fitted CHO, `z_det` = 3.0: fraction of planned lesions at or
above `z_det` (median `z_ref`), thick slices as acquired.

| sequence, contrast | 10 mm^3 | 27 mm^3 | 64 mm^3 | 100 mm^3 | 200 mm^3 |
|---|---|---|---|---|---|
| AXFLAIR +0.24 | 0 % (z 0.3) | 3 % (z 0.4) | 5 % (z 0.5) | 3 % (z 0.5) | 3 % (z 0.7) |
| AXFLAIR +1.00 | 15 % (z 1.6) | 26 % (z 2.1) | 35 % (z 2.6) | 46 % (z 2.9) | 63 % (z 3.4) |
| AXFLAIR +1.50 | 38 % (z 2.6) | 57 % (z 3.2) | 76 % (z 4.0) | 79 % (z 4.5) | 85 % (z 5.1) |
| AXT1 -0.21 | 2 % (z 0.6) | 1 % (z 0.8) | 1 % (z 0.8) | 0 % (z 0.7) | 0 % (z 0.8) |
| AXT1 -0.60 | 4 % (z 1.4) | 5 % (z 1.7) | 5 % (z 1.8) | 3 % (z 1.6) | 7 % (z 1.7) |
| AXT1 -1.00 | 27 % (z 2.3) | 44 % (z 2.8) | 41 % (z 2.8) | 30 % (z 2.6) | 44 % (z 2.8) |

Expected eligible fraction: 18 % of planned insertions measured on the full run (the patch-only check estimated 25 %; the full run cross-fits per slice unit and is stricter, so 18 % is the number to size with) (FLAIR
36 %, T1 14 %), against 3 % for the pre-registered bank.
T1 does not reach a majority at any level, because -1.0 is already a signal
void; this is a property of the thick-slice data and is stated as such.

### Option A1a: a leaner confirmatory grid (for the authors to decide)

The subtle strata (FLAIR 0.24, T1 -0.21) and the 10 mm^3 volume are almost
entirely sub-threshold. Keeping every combination spends about 40 % of the
GPU reconstruction budget on pairs that cannot enter the erasure endpoint.
A1a keeps one sub-threshold contrast control (the subtle level at 100 mm^3
only) and one sub-threshold volume control (10 mm^3 at the top level only),
and runs the full contrast axis at 27, 64, 100 and 200 mm^3. The
model-independent characterisation already computed for the full grid is
unaffected either way.

## Amendment A2: H1's comparator (an endpoint's test, not the endpoint)

The classical audit at R = 8 on the amended bank found erasure rates of 0.68
(zero-filling) and 0.87 (CG-SENSE at the audit's noise-set lambda) among
eligible pairs, with no learned prior involved (`docs/VALIDATION.md`,
2026-09-11). The pre-registered H1 — "the silent erasure rate exceeds
`pi0` = 0.05 in at least one model family" — is therefore satisfied by the
acquisition alone and carries no information about learned reconstruction.

**The endpoint does not change.** It remains the silent erasure rate, defined
exactly as pre-registered. What changes is the comparator:

| | Pre-registered | Amended (A2) |
| --- | --- | --- |
| Endpoint | silent erasure rate | unchanged |
| Test | one-sided against `pi0` = 0.05, per model | two-sided paired contrast against a no-prior baseline, per model |
| Estimand | the rate | `Delta` = rate(model) - rate(baseline), same lesions, acquisitions and seeds |
| `pi0` | the test | reported as a descriptive reference for clinical relevance |
| Strata | volume, contrast, R | unchanged, plus the measured fraction `mu_lambda` |
| Sizing | half-width 0.02 on a rate | half-width 0.05 on `Delta`; not more demanding, because a paired contrast of correlated rates has smaller variance than either rate |

The comparator is the fourth arm added to Experiment 1: a Tikhonov-regularised
sensitivity-encoded reconstruction at an image-quality-oriented regularisation
level (`cg_sense_tuned`). Pairing removes everything the two arms share, so the
contrast isolates the prior's contribution, which a threshold test cannot do
once the acquisition is known to erase lesions on its own.

Reporting at matched `mu_lambda` is part of the amendment: a lesion the
acquisition barely records cannot be preserved by any method. R = 4, where
`mu_lambda` is essentially one, is the condition in which a difference is most
cleanly attributable to the prior, and the contrast is reported by R as well
as pooled.

## Amendment A2b: reporting strata (a reporting addition, not a test change)

The classical audit across accelerations (`docs/VALIDATION.md`, 2026-09-12)
showed that absolute erasure rates are governed by how far a lesion's
reference detectability sits above `z_det`, because the pre-registered
`z_det` = 3.0 and `z_miss` = 2.0 are a factor 1.5 apart while undersampling
alone costs a factor 2.0 to 2.8 in detectability. Erasure rates for lesions in
the band `z_ref` 3.0 to 3.5 ran at 66 to 89 %, against 24 to 41 % for `z_ref`
above 6, under classical reconstruction with no prior at all.

The paired contrast of A2 is unaffected, because the loss is common to both
arms and cancels. What is added here is reporting: the contrast is reported
stratified by `z_ref` band (3.0-3.5, 3.5-4.5, 4.5-6, above 6) as well as by
`mu_lambda` and acceleration. No threshold, definition or test changes.

A reference line is also fixed in advance: the tuned classical arm tracks the
SNR-limited prediction `1/sqrt(R)` for detectability loss, so the confirmatory
result is read as whether a learned prior recovers detectability beyond the
SNR limit, rather than against an arbitrary rate.

## Amendment A3: scope of the claims

No public raw-k-space benchmark contains acute infarction or haemorrhage, and
the focal lesions annotated in fastMRI+ are predominantly chronic: the 113
boxes labelled lacunar infarct have a measured core contrast of median
magnitude near 0.2, against the bank's conspicuous levels of about four times
that. Those levels are an assumption about how an acute lesion would appear,
not a measurement of one. All claims are therefore stated as being about
focal structure of specified size and contrast under a specified acquisition.
The clinical findings that motivate the work are named as motivation only.
This is a scope statement rather than a change of procedure, recorded here so
that the limit is on the record before any confirmatory result exists.

## Realised sample size (to be filled after the pilot)

The pilot estimates `rho` and nothing else; its subjects are excluded from
the confirmatory set. Pairs per model-by-R cell follow the pre-registered
rule with `pi0` = 0.05, half-width 0.02 and the design effect
`1 + (m - 1) rho`. Note from the statistics simulation: H2's power depends on
the number of *erased* pairs, so it is tested on the pooled confirmatory
sample, and if the realised erasure rate is far below `pi0` the addendum
should say H2 is under-powered before the run.

## Provenance

Code at the time of this draft: `csa-audit` commit (fill), Zenodo concept DOI
10.5281/zenodo.22648224. Data: fastMRI brain multi-coil val batches 0 to 2
and train batch 0 (train volumes used for contrast statistics only, never as
hosts); fastMRI+ brain annotations.
