#!/usr/bin/env python3
"""Physics demonstration behind Experiment 2, with a classical reconstructor.

VALIDATION, NOT A RESULT. Experiment 2 proper decomposes the errors of the
learned models' erased pairs, which need a GPU. What can be shown now, on real
acquisitions, is the operator-level content of Proposition 1:

  (a) the change in the k-space residual produced by a pure null-space edit
      delta = Q_lambda(ell') is bounded by ||A delta|| and sits far below the
      residual itself;
  (b) removing the null-space component Q(e) of a reconstruction's error
      leaves the residual essentially unchanged, i.e. the residual tracks
      P(e) only;
  (c) the lesion transfer of CG-SENSE has t_R near 1 and t_N near 0.

Rows come from the bank characterization output so the same lesions are used.
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
from csa.nullspace.projectors import null_filter, range_filter, transfer_statistics  # noqa: E402
from csa.physics.masks import default_center_fraction, equispaced_mask  # noqa: E402
from csa.physics.operator import SenseOperator, ifft2c, sense_combine  # noqa: E402
from csa.recon.base import CGSense  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--characterization", default="outputs/bank_characterization.csv")
    ap.add_argument("--data", required=True)
    ap.add_argument("--maps-cache", default="/Volumes/T9/fastMRI_brain/cache/maps")
    ap.add_argument("--out", default="outputs/exp2_physics_demo.csv")
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--acceleration", type=int, default=8)
    ap.add_argument("--cg-tol", type=float, default=1e-5)
    ap.add_argument("--cg-maxiter", type=int, default=80)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    ch = pd.read_csv(args.characterization).rename(columns={"lambda": "lam"})  # "lambda" is a keyword; itertuples renames it
    ch = ch[(ch["acceleration"] == args.acceleration) & (ch["mask_type"] == "equispaced")]
    ch = ch.drop_duplicates(["file", "slice", "site_row", "site_col", "volume_mm3"])
    rng = np.random.default_rng(args.seed)
    picks = ch.sample(n=min(args.n, len(ch)), random_state=args.seed)
    files = {os.path.splitext(f)[0]: os.path.join(dp, f) for dp, _, fs in os.walk(args.data) for f in fs if f.endswith(".h5")}
    rows = []
    for i, r in enumerate(picks.itertuples(), 1):
        t0 = time.time()
        path = files[r.file]
        with h5py.File(path, "r") as f:
            ks = np.asarray(f["kspace"][int(r.slice)]).astype(np.complex128)
        maps = np.load(os.path.join(args.maps_cache, f"{r.file}_s{int(r.slice):02d}.npz"))["maps"].astype(np.complex128)
        X = sense_combine(ifft2c(ks), maps)
        mask = equispaced_mask(ks.shape[-1], args.acceleration, default_center_fraction(args.acceleration))
        op = SenseOperator(maps, mask)
        lam = float(r.lam)
        contrast = float(ch[(ch["file"] == r.file) & (ch["slice"] == r.slice)]["contrast"].median())
        ell = lesion_perturbation(X, (int(r.site_row), int(r.site_col)), radius_px=float(r.radius_px), contrast=contrast, fill=float(r.fill))
        y = op.undersample(ks); y_cf = y + op.forward(ell)
        rec = CGSense(lam=lam, tol=args.cg_tol, maxiter=args.cg_maxiter)
        Xh = rec(y, mask, maps); Xh_cf = rec(y_cf, mask, maps)
        X_cf = X + ell
        e = Xh_cf - X_cf
        Pe = range_filter(op, e, lam, tol=args.cg_tol, maxiter=args.cg_maxiter); Qe = e - Pe
        ny = np.linalg.norm(y_cf)
        r_full = np.linalg.norm(op.forward(Xh_cf) - y_cf) / ny
        r_minus_null = np.linalg.norm(op.forward(Xh_cf - Qe) - y_cf) / ny
        r_minus_range = np.linalg.norm(op.forward(Xh_cf - Pe) - y_cf) / ny
        # a pure null-space edit from a second lesion elsewhere
        a = np.abs(X); ssq = (np.abs(maps) ** 2).sum(0)
        cand = np.argwhere((ssq > 0.5) & (a > 0.25 * a[ssq > 0.5].max()))
        s2 = tuple(int(v) for v in cand[rng.integers(len(cand))])
        ell2 = lesion_perturbation(X, s2, radius_px=float(r.radius_px), contrast=contrast, fill=float(r.fill))
        delta = null_filter(op, ell2, lam, tol=args.cg_tol, maxiter=args.cg_maxiter)
        r_edit = np.linalg.norm(op.forward(Xh_cf + delta) - y_cf) / ny
        tR, tN = transfer_statistics(op, Xh_cf - Xh, ell, lam, tol=args.cg_tol, maxiter=args.cg_maxiter)
        rows.append({"file": r.file, "slice": int(r.slice), "acceleration": args.acceleration, "volume_mm3": float(r.volume_mm3),
                     "contrast": contrast, "lambda": lam,
                     "err_energy_total": float(np.vdot(e, e).real), "err_energy_range": float(np.vdot(Pe, Pe).real),
                     "err_energy_null": float(np.vdot(Qe, Qe).real),
                     "null_error_share": float(np.vdot(Qe, Qe).real / np.vdot(e, e).real),
                     "residual_full": r_full, "residual_without_null_error": r_minus_null,
                     "residual_without_range_error": r_minus_range,
                     "edit_norm_rel": float(np.linalg.norm(delta) / np.linalg.norm(Xh_cf)),
                     "edit_A_norm_rel": float(np.linalg.norm(op.forward(delta)) / ny),
                     "residual_with_edit": r_edit, "residual_change_from_edit": abs(r_edit - r_full),
                     "bound_holds": bool(abs(r_edit - r_full) <= np.linalg.norm(op.forward(delta)) / ny + 1e-12),
                     "edit_invisibility": float(np.linalg.norm(op.forward(delta)) / np.linalg.norm(delta)),
                     "t_R_cgsense": tR, "t_N_cgsense": tN, "seconds": round(time.time() - t0, 1)})
        print(f"[{i}/{len(picks)}] {r.file} s{int(r.slice)} vol {r.volume_mm3}: residual {r_full:.4f} -> "
              f"minus null err {r_minus_null:.4f}, minus range err {r_minus_range:.4f}; null-edit changes residual by "
              f"{abs(r_edit-r_full):.2e} (bound {np.linalg.norm(op.forward(delta))/ny:.2e}); t_R {tR:.3f} t_N {tN:.3f}", flush=True)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    print("\nSUMMARY (classical reconstructor; validation, not Experiment 2's result)")
    print(df[["null_error_share", "residual_full", "residual_without_null_error", "residual_without_range_error",
              "residual_change_from_edit", "edit_A_norm_rel", "edit_invisibility", "t_R_cgsense", "t_N_cgsense"]].describe().T[["mean", "min", "max"]])
    print("Proposition 1 bound held in", int(df["bound_holds"].sum()), "of", len(df), "cases")


if __name__ == "__main__":
    main()
