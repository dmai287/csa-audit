#!/usr/bin/env python3
"""Measurement-side part of Experiment 3 (acquisition-shift stress), no model.

Experiment 3 proper compares the stability of observer detectability and t_N
against PSNR/SSIM under acquisition shift; that needs the learned models. The
model-independent part is how much of a lesion the acquisition measures
(kappa, mu_lambda) as mask pattern, noise level and coil configuration are
intervened on at fixed anatomy and lesion. Runs on CPU from the bank.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.lesion.insert import lesion_perturbation  # noqa: E402
from csa.nullspace.projectors import cg_solve, measurement_gain, real_inner  # noqa: E402
from csa.physics.masks import default_center_fraction, equispaced_mask, random_mask  # noqa: E402
from csa.physics.operator import SenseOperator, ifft2c, sense_combine  # noqa: E402


def mu_of(op, ell, lam, tol, maxiter):
    z, info = cg_solve(lambda v: op.normal(v) + lam * v, op.normal(ell), tol=tol, maxiter=maxiter)
    return real_inner(z, z) / real_inner(ell, ell), info["iterations"], info["relative_residual"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--characterization", default="outputs/bank_characterization.csv")
    ap.add_argument("--data", required=True)
    ap.add_argument("--maps-cache", default="/Volumes/T9/fastMRI_brain/cache/maps")
    ap.add_argument("--out", default="outputs/exp3_measurement_side.csv")
    ap.add_argument("--n-lesions", type=int, default=12)
    ap.add_argument("--acceleration", type=int, default=8)
    ap.add_argument("--cg-tol", type=float, default=1e-5)
    ap.add_argument("--cg-maxiter", type=int, default=80)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    ch = pd.read_csv(args.characterization).rename(columns={"lambda": "lam"})  # "lambda" is a keyword; itertuples renames it
    ch = ch[(ch["acceleration"] == args.acceleration) & (ch["mask_type"] == "equispaced")]
    ch = ch.drop_duplicates(["file", "slice", "site_row", "site_col", "volume_mm3"])
    picks = ch.sample(n=min(args.n_lesions, len(ch)), random_state=args.seed)
    files = {os.path.splitext(f)[0]: os.path.join(dp, f) for dp, _, fs in os.walk(args.data) for f in fs if f.endswith(".h5")}
    R = args.acceleration; cf = default_center_fraction(R)
    rows = []
    for i, r in enumerate(picks.itertuples(), 1):
        with h5py.File(files[r.file], "r") as f:
            ks = np.asarray(f["kspace"][int(r.slice)]).astype(np.complex128)
        maps = np.load(os.path.join(args.maps_cache, f"{r.file}_s{int(r.slice):02d}.npz"))["maps"].astype(np.complex128)
        X = sense_combine(ifft2c(ks), maps)
        ssq = (np.abs(maps) ** 2).sum(0); support = ssq > 0.5
        lam0 = float(r.lam)
        ell = lesion_perturbation(X, (int(r.site_row), int(r.site_col)), radius_px=float(r.radius_px), contrast=1.0, fill=float(r.fill))
        n_cols, n_coils = ks.shape[-1], ks.shape[0]
        conditions = []
        for off in range(min(R, 4)):
            conditions.append((f"equispaced_offset_{off}", equispaced_mask(n_cols, R, cf, offset=off), maps, lam0))
        for s in range(3):
            conditions.append((f"random_seed_{s}", random_mask(n_cols, R, cf, seed=s), maps, lam0))
        for k in (2.0, 4.0):
            conditions.append((f"noise_x{k:g}", equispaced_mask(n_cols, R, cf), maps, lam0 * k * k))
        if n_coils >= 4:
            for name, idx in (("coils_first_half", np.arange(n_coils // 2)), ("coils_second_half", np.arange(n_coils // 2, n_coils)),
                              ("coils_alternate", np.arange(0, n_coils, 2))):
                sub = maps[idx]
                norm = np.sqrt((np.abs(sub) ** 2).sum(0)); norm[norm == 0] = 1.0
                conditions.append((name, equispaced_mask(n_cols, R, cf), sub / norm, lam0))
        for name, mask, m_used, lam in conditions:
            t0 = time.time()
            op = SenseOperator(m_used, mask)
            kappa = measurement_gain(op, ell)
            mu, it, rr = mu_of(op, ell, lam, args.cg_tol, args.cg_maxiter)
            rows.append({"file": r.file, "slice": int(r.slice), "site_row": int(r.site_row), "site_col": int(r.site_col),
                         "volume_mm3": float(r.volume_mm3), "acceleration": R, "condition": name,
                         "n_coils_used": int(m_used.shape[0]), "sampled_fraction": float(mask.mean()), "lambda": lam,
                         "kappa": kappa, "kappa2": kappa ** 2, "mu_lambda": mu, "cg_iterations": it, "cg_rel_residual": rr,
                         "seconds": round(time.time() - t0, 1)})
        print(f"[{i}/{len(picks)}] {r.file} s{int(r.slice)} vol {r.volume_mm3}: {len(conditions)} conditions done", flush=True)
        pd.DataFrame(rows).to_csv(args.out, index=False)
    df = pd.DataFrame(rows)
    print("\nSUMMARY: mu_lambda by condition (mean over lesions)")
    print(df.groupby("condition")[["kappa2", "mu_lambda", "sampled_fraction"]].mean().round(4))


if __name__ == "__main__":
    main()
