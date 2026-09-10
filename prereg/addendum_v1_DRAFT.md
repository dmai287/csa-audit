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

Expected eligible fraction: 25 % of planned insertions (FLAIR
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
