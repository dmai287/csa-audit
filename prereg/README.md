# Pre-registration

The statistical analysis plan (paper Section 4.8) and the thresholds in
`configs/thresholds.yaml` are deposited on OSF as a time-stamped
pre-registration **before any reconstruction is run**. Any later deviation is
logged in `CHANGELOG.md` with its date and reason.

## v1

- Source: `prereg_v1.tex` (compile with pdflatex; TeX Live)
- Deposited PDF: `prereg_v1.pdf` (5 pages)
- SHA-256: `720872a701ff16fe5de2059525d3c24a6cdc2a155b58a204f61e4b9769cc5a03`
  (also in `prereg_v1.sha256`; verify with `shasum -a 256 -c prereg_v1.sha256`)
- Prepared: 2026-09-08
- OSF identifier: (pending deposit)
- Deposit date: (pending deposit)

Status at preparation: no reconstruction had been run, no pilot had been run,
and no experimental result of any kind existed. The thresholds were fixed from
convention and prior literature with no contact with reconstruction output.

The code state matching this document is Zenodo concept DOI
`10.5281/zenodo.22648224` (release `v0.1.0`, version DOI
`10.5281/zenodo.22648225`).

## Thresholds fixed by v1

| Symbol | Value | Meaning |
| --- | --- | --- |
| `z_det` | 3.0 | Reference detectability required for a lesion to enter the erasure analysis |
| `z_miss` | 2.0 | Reconstruction detectability below which the lesion counts as missed |
| `z_fs` | 4.0 | False-structure threshold, in sd of null-space background error |
| `pi0` | 0.05 | Silent erasure rate taken as the floor of clinical relevance (H1) |

The gap between `z_miss` and `z_det` is deliberate: pairs falling in `[2.0,
3.0)` are indeterminate and count as neither erased nor preserved, so
borderline pairs cannot inflate the erasure rate.

## Outstanding: the sizing addendum

v1 fixes the sample-size *rule* but not the realised N, because the
intra-subject correlation `rho` is unknown until the 20-subject design pilot
runs. v1 commits to four constraints on that pilot: it estimates `rho` and
nothing else, no endpoint is computed on it, its subjects are excluded from the
confirmatory set, and the realised pairs-per-cell is posted as a dated OSF
addendum **before** any confirmatory reconstruction.

So the order is: deposit v1 → run the pilot → post the addendum → run the
confirmatory experiments.
