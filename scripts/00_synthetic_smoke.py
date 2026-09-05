#!/usr/bin/env python3
"""End-to-end smoke run of the audit on a synthetic phantom (no fastMRI data).

Builds a phantom scene, fits the observer on lesion-absent/present training
pairs, then audits inserted lesions under zero-filled and CG-SENSE
reconstruction at several accelerations. Output: one CSV row per pair. This
exercises every module the real experiments will use; its numbers mean
nothing scientifically.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.audit.run import PairSpec, audit_pair  # noqa: E402
from csa.lesion.insert import lesion_perturbation  # noqa: E402
from csa.observer.cho import CHO, extract_roi, laguerre_gauss_channels  # noqa: E402
from csa.physics.masks import default_center_fraction, equispaced_mask  # noqa: E402
from csa.physics.operator import SenseOperator  # noqa: E402
from csa.recon.base import CGSense, ZeroFilled  # noqa: E402
from csa.synthetic import coil_maps, multicoil_kspace, phantom, support_mask  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=96)
    ap.add_argument("--pairs", type=int, default=6)
    ap.add_argument("--noise", type=float, default=0.02)
    ap.add_argument("--out", default="outputs/smoke.csv")
    ap.add_argument("--roi", type=int, default=24)
    args = ap.parse_args()

    rng = np.random.default_rng(0)
    X = phantom(args.n)
    maps = coil_maps(args.n, 6)
    y_full = multicoil_kspace(X, maps, noise_sigma=args.noise, seed=1)
    support = support_mask(X)
    rows_idx = np.argwhere(support)
    lam = args.noise ** 2 / np.var(np.abs(X))
    U = laguerre_gauss_channels(args.roi, 6, a_px=args.roi / 3)
    recons = {"zero_filled": ZeroFilled(), "cg_sense": CGSense(lam=lam, tol=1e-8)}
    rows = []
    for R in (4, 8):
        mask = equispaced_mask(args.n, R, default_center_fraction(R))
        op = SenseOperator(maps, mask)
        y = op.undersample(y_full)
        for name, rec in recons.items():
            # lesion-absent ensemble for the observer: reconstructions of the factual arm at random sites
            Xh0 = rec(y, mask, maps)
            sites = rows_idx[rng.choice(len(rows_idx), size=40, replace=False)]
            absent = np.stack([extract_roi(Xh0, tuple(s), args.roi) for s in sites])
            present = []
            for s in sites:
                ell = lesion_perturbation(X, tuple(s), radius_px=2.5, contrast=0.6)
                present.append(extract_roi(rec(y + op.forward(ell), mask, maps), tuple(s), args.roi))
            present = np.stack(present)
            cho = CHO(U).fit(present[:20], absent[:20])
            t0 = cho.decide(absent[20:])
            # reference-domain observer: same channels, template fitted on reference patches
            absent_ref = np.stack([extract_roi(X, tuple(s), args.roi) for s in sites])
            present_ref = np.stack([extract_roi(X + lesion_perturbation(X, tuple(s), 2.5, 0.6), tuple(s), args.roi)
                                    for s in sites])
            cho_ref = CHO(U).fit(present_ref[:20], absent_ref[:20])
            t0_ref = cho_ref.decide(absent_ref[20:])
            for k in range(args.pairs):
                s = tuple(rows_idx[rng.integers(len(rows_idx))])
                contrast = float(rng.choice([0.3, 0.6, 1.0]))
                ell = lesion_perturbation(X, s, radius_px=2.5, contrast=contrast)
                spec = PairSpec("synthetic", 0, int(s[0]), int(s[1]), 27.0, contrast, name, R, "equispaced", 0)
                rows.append(audit_pair(spec, op, y_full, X, ell, rec, lam, args.roi, cho, t0, t0_ref,
                                       radius_px=2.5, cho_ref=cho_ref))
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    cols = ["model", "acceleration", "contrast", "mu_lambda", "t_R", "t_N", "z_ref", "z_recon", "psnr_cf", "residual_cf"]
    print(df[cols].round(3).to_string(index=False))
    print(f"wrote {len(df)} rows to {args.out}")


if __name__ == "__main__":
    main()
