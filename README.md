# csa-audit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22648224.svg)](https://doi.org/10.5281/zenodo.22648224)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Code for **The Causal Safety Audit: Counterfactual Stress-Testing for
Deep-Learning MRI Reconstruction** (Mai, Pham, Nguyen; manuscript in
preparation, target *Physics in Medicine & Biology*). Companion to the
systematic review *The Metric Blind Spot in Deep Learning-Based Brain MRI
Reconstruction* (SSRN preprint, doi:10.2139/ssrn.7362258).

> **Status (2026-09-05): no experiment has been run.** This repository holds
> the audit implementation and its tests. Every number the paper will report
> is to be produced by the scripts here from a pre-registered plan; none exists
> yet.

## What the audit does

For a lesion-free acquisition `Y = A X + eta`, a lesion perturbation `ell` is
inserted on the measurement side of the forward operator, `Y_cf = Y + A ell`,
so the factual and counterfactual arms share the measured noise realisation,
sampling mask and coil sensitivities exactly. Both arms are reconstructed with
the same model and seed, and the pair is read out on three instruments:

1. **Null-space statistics.** Filtered projectors `P = (A^H A + lam I)^-1 A^H A`
   and `Q = I - P` split the lesion transfer `dX = X_cf_hat - X_hat` into the
   part the measurement constrained and the part the prior supplied: the lesion
   measurement fraction `mu`, the measured transfer `t_R = <A dX, A ell>/||A ell||^2`
   (exactly 1 under hard data consistency, which is all a residual gate checks)
   and the null-space transfer `t_N` on `Q`.
2. **Task-based detectability.** A channelised Hotelling observer with
   Laguerre-Gauss channels, ensemble `d'` and AUC, and a per-lesion z-score.
3. **Image scores.** PSNR and SSIM in the fastMRI convention, ROI CNR, ROI SSIM.

Erasure and *silent* erasure (lesion lost while global PSNR stays inside its
lesion-free range) are defined before any run from `configs/thresholds.yaml`.

## Layout

```
csa/physics     SENSE forward operator and adjoint, masks, ESPIRiT wrapper
csa/lesion      Algorithm 1: additive insertion carried through A
csa/nullspace   CG filtered projectors, mu / t_R / t_N, false-structure count
csa/observer    channelised Hotelling observer
csa/metrics     PSNR, SSIM, ROI measures
csa/audit       one pair -> one row; erasure definitions
csa/recon       zero-filled, CG-SENSE; adapters for VarNet, complex U-Net, diffusion
csa/stats       cluster bootstrap, TOST equivalence, AUC intervals, pairs-per-cell rule
csa/io          fastMRI HDF5 access and integrity report
scripts/        smoke run, integrity check, bank, pilot, experiments, statistics
tests/          gates on synthetic data (no fastMRI files needed)
configs/        masks, lesion bank, thresholds (null until pre-registered), models
prereg/         the deposited pre-registration PDF and its hash
docs/           model notes and the steps that need a browser
```

## Install and test

```
pip install -e ".[test]"          # NumPy core, works on Python 3.8+
pip install -e ".[recon]"         # torch, fastmri, sigpy for the learned models
python -m pytest -q
python scripts/00_synthetic_smoke.py --out outputs/smoke.csv
python scripts/08_stats.py --table outputs/smoke.csv
```

The tests are the gates from the build plan: the adjoint identity of the
operator; equality of the two lesion-insertion paths to machine precision;
closed-form behaviour of the filtered projectors on a single-coil operator,
their Hermitian symmetry, the transfer identities, and the fact that a linear
data-consistent reconstructor transfers exactly `P ell`; the residual bound of
Proposition 1 on a real operator; and that the observer never beats the ideal
linear observer and scales linearly with contrast.

## Workflow

1. `scripts/01_integrity.py` on the fastMRI brain multi-coil files: volumes with
   zero-filled k-space are excluded, spacing and coil counts recorded.
2. `scripts/02_bank.py` builds the lesion-bank manifest from fastMRI+.
3. Deposit the pre-registration (see `prereg/`), copy thresholds into
   `configs/thresholds.yaml`, never edit them afterwards.
4. `scripts/03_pilot.py` on a 20-subject pilot fixes pairs per cell.
5. `scripts/04_exp1.py` onwards produce the per-pair table with provenance
   columns (git commit, config hash, checkpoint hash, data checksums).
6. `scripts/08_stats.py` and `scripts/09_figures.py` read only that table.

Steps that need a browser (Zenodo integration, fastMRI data agreement, OSF,
ethics) are listed in `docs/BROWSER_STEPS.md`. Model notes are in
`docs/MODELS.md`. Deviations from the pre-registration go in `CHANGELOG.md`.

## Rules

- Models are audited as released; nothing is retrained on the lesion bank.
- No lesion-bearing slice enters any training or calibration step.
- Both arms of a pair share the seed; seeds are then varied across replicates.
- Data and weights are never committed (`data/`, `checkpoints/`, `outputs/`
  are ignored).

## License and citation

Code under the MIT License. Cite via `CITATION.cff`. Data are fastMRI and
fastMRI+, under their own terms.

Releases are archived on Zenodo. Cite the concept DOI
[`10.5281/zenodo.22648224`](https://doi.org/10.5281/zenodo.22648224), which
always resolves to the most recent version; cite a version DOI (`v0.1.0` is
[`10.5281/zenodo.22648225`](https://doi.org/10.5281/zenodo.22648225)) when a
specific release must be pinned.
