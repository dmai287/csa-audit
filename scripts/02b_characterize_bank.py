#!/usr/bin/env python3
"""Model-independent half of Experiment 1, computable without any reconstructor.

For every planned insertion in the bank manifest this computes, from the fully
sampled reference and the real forward operator only:

  kappa(ell)      measurement gain ||A ell|| / ||ell||          (per site, volume, R)
  mu_lambda(ell)  measured fraction ||P_lambda ell||^2/||ell||^2  (per site, volume, R)
  z_ref           per-lesion detectability of the lesion in the fully sampled
                  reference, from a channelised Hotelling observer fitted on
                  reference patches with two-fold cross-fitting by file
  roi_cnr_ref     ROI contrast-to-noise of the inserted lesion in the reference

No reconstruction output of any model under audit is touched, no erasure
indicator is computed, and the pre-registered thresholds are not applied.
These quantities are needed by the confirmatory analysis regardless, and they
depend only on the lesion and the acquisition, so computing them once here
saves recomputing them for each of the three models later.

Pass 1 (multiprocess, one unit per (file, slice)): ESPIRiT maps (cached),
reference image, noise level, and for each site/volume/R the CG solves.
Reference patches are saved so the observer can be fitted in pass 2 across
files. Resumable: (file, slice) units already in the output are skipped.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import h5py
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.lesion.insert import lesion_perturbation, slice_fill_fraction, soft_disc, sphere_radius_mm  # noqa: E402
from csa.metrics.image import annulus, roi_cnr  # noqa: E402
from csa.nullspace.projectors import cg_solve, measurement_gain, real_inner  # noqa: E402
from csa.observer.cho import CHO, extract_roi, laguerre_gauss_channels, per_lesion_z  # noqa: E402
from csa.physics.espirit import espirit_maps  # noqa: E402
from csa.physics.masks import default_center_fraction, equispaced_mask, random_mask  # noqa: E402
from csa.physics.operator import SenseOperator, ifft2c, sense_combine  # noqa: E402


def load_maps(path, sl, cache_dir, calib_width=24):
    key = f"{os.path.splitext(os.path.basename(path))[0]}_s{sl:02d}.npz"
    cp = os.path.join(cache_dir, key)
    if os.path.exists(cp):
        return np.load(cp)["maps"].astype(np.complex128)
    with h5py.File(path, "r") as f:
        ks = np.asarray(f["kspace"][sl])
    maps = espirit_maps(ks, calib_width=calib_width)
    os.makedirs(cache_dir, exist_ok=True)
    np.savez_compressed(cp, maps=maps.astype(np.complex64))
    return maps.astype(np.complex128)


def measured_fraction_with_info(op, ell, lam, tol, maxiter):
    b = op.normal(ell)
    z, info = cg_solve(lambda v: op.normal(v) + lam * v, b, tol=tol, maxiter=maxiter)
    mu = real_inner(z, z) / real_inner(ell, ell)
    return mu, info


def process_unit(job):
    """One (file, slice): returns (rows, patches)."""
    (path, sl, rows_spec, args) = job
    t_unit = time.time()
    with h5py.File(path, "r") as f:
        ks = np.asarray(f["kspace"][sl]).astype(np.complex128)
    maps = load_maps(path, sl, args["maps_cache"])
    X = sense_combine(ifft2c(ks), maps)
    ssq = (np.abs(maps) ** 2).sum(0)
    support = ssq > 0.5
    corner = ks[:, :16, :16]
    sigma_k = float(np.sqrt(0.5 * (corner.real.var() + corner.imag.var())))
    lam = float(sigma_k ** 2 / np.var(np.abs(X)[support]))
    stem = os.path.splitext(os.path.basename(path))[0]
    spec0 = rows_spec[0]
    sx = float(spec0.get("spacing_x_mm") or 0.7)
    thick = float(spec0.get("thickness_mm") or 5.0)
    roi = args["roi"]
    rng = np.random.default_rng(hash((stem, sl)) % (2 ** 32))

    ops = {}
    for R in args["accelerations"]:
        for mt in args["mask_types"]:
            cf = default_center_fraction(R)
            m = equispaced_mask(ks.shape[-1], R, cf) if mt == "equispaced" else random_mask(ks.shape[-1], R, cf, seed=0)
            ops[(R, mt)] = SenseOperator(maps, m)

    out_rows, patches = [], []
    # unique (site, volume) combinations; contrasts share them
    seen = {}
    for spec in rows_spec:
        key = (int(spec["site_row"]), int(spec["site_col"]), float(spec["volume_mm3"]))
        seen.setdefault(key, []).append(float(spec["contrast"]))
    for (row, col, vol), contrasts in seen.items():
        rmm = sphere_radius_mm(vol); rpx = rmm / sx; fill = slice_fill_fraction(rmm, thick)
        site = (row, col)
        try:
            ell_unit = lesion_perturbation(X, site, radius_px=rpx, contrast=1.0, fill=fill)
        except ValueError:
            continue
        e_unit = float(real_inner(ell_unit, ell_unit))
        if e_unit == 0.0:
            continue
        disc = soft_disc(X.shape, site, rpx, taper_px=0.0) > 0.5
        ring = annulus(X.shape, site, rpx + 2, rpx + 8) & support
        absent_patch = extract_roi(X, site, roi)
        ell_patch = extract_roi(ell_unit, site, roi)  # |ell| patch; present = absent + c * ell_patch (phase-coherent)
        patches.append({"file": stem, "slice": sl, "sequence": spec0["sequence"], "site_row": row, "site_col": col,
                        "volume_mm3": vol, "radius_px": rpx, "absent": absent_patch.astype(np.float32),
                        "ell": ell_patch.astype(np.float32)})
        for (R, mt), op in ops.items():
            kappa = measurement_gain(op, ell_unit)
            t0 = time.time()
            mu, info = measured_fraction_with_info(op, ell_unit, lam, args["cg_tol"], args["cg_maxiter"])
            base = {"file": stem, "slice": sl, "sequence": spec0["sequence"], "site_row": row, "site_col": col,
                    "volume_mm3": vol, "radius_mm": rmm, "radius_px": rpx, "fill": fill,
                    "acceleration": R, "mask_type": mt, "sampled_fraction": float(op.mask.mean()),
                    "n_coils": int(ks.shape[0]), "noise_sigma_k": sigma_k, "lambda": lam,
                    "kappa": kappa, "kappa2": kappa ** 2, "mu_lambda": mu,
                    "cg_iterations": info["iterations"], "cg_rel_residual": info["relative_residual"],
                    "cg_seconds": round(time.time() - t0, 1), "lesion_energy_unit": e_unit}
            for c in contrasts:
                Xc = X + c * ell_unit
                out_rows.append({**base, "contrast": c,
                                 "roi_cnr_ref": roi_cnr(Xc, disc, ring) if ring.any() else np.nan,
                                 "lesion_energy": e_unit * c * c})
    # extra lesion-absent reference patches for the observer's covariance estimate
    a = np.abs(X); cand = np.argwhere(support & (a > 0.25 * a[support].max()))
    for k in range(args["extra_absent"]):
        r, c = cand[rng.integers(len(cand))]
        patches.append({"file": stem, "slice": sl, "sequence": spec0["sequence"], "site_row": int(r), "site_col": int(c),
                        "volume_mm3": -1.0, "radius_px": 0.0, "absent": extract_roi(X, (r, c), roi).astype(np.float32),
                        "ell": None})
    for rr in out_rows:
        rr["unit_seconds"] = round(time.time() - t_unit, 1)
    return out_rows, patches


def fit_reference_observer(patches, rows, roi, seed=0):
    """Two-fold cross-fitted CHO per (sequence, volume): z_ref for every planned lesion."""
    z = {}
    dprime_rows = []
    rng = np.random.default_rng(seed)
    df_p = pd.DataFrame([{k: v for k, v in p.items() if k not in ("absent", "ell")} for p in patches])
    for seq in sorted(df_p["sequence"].unique()):
        files = sorted(df_p.loc[df_p["sequence"] == seq, "file"].unique())
        if len(files) < 2:
            continue
        rng.shuffle(files)
        folds = [set(files[::2]), set(files[1::2])]
        absent_all = [p for p in patches if p["sequence"] == seq]
        for vol in sorted(v for v in df_p.loc[df_p["sequence"] == seq, "volume_mm3"].unique() if v > 0):
            lesion_p = [p for p in patches if p["sequence"] == seq and p["volume_mm3"] == vol]
            if not lesion_p:
                continue
            rpx = float(np.median([p["radius_px"] for p in lesion_p]))
            U = laguerre_gauss_channels(roi, 6, a_px=float(np.clip(3.0 * rpx, 4.0, roi / 2)))
            contrasts = sorted({r["contrast"] for r in rows if r["sequence"] == seq and r["volume_mm3"] == vol})
            c_fit = contrasts[len(contrasts) // 2]
            for k in (0, 1):
                train, test = folds[k], folds[1 - k]
                tr_abs = np.stack([p["absent"] for p in absent_all if p["file"] in train])
                tr_pres = np.stack([p["absent"] + c_fit * p["ell"] for p in lesion_p if p["file"] in train])
                te_abs = np.stack([p["absent"] for p in absent_all if p["file"] in test])
                if len(tr_abs) < 8 or len(tr_pres) < 4 or len(te_abs) < 8:
                    continue
                cho = CHO(U).fit(tr_pres.reshape(len(tr_pres), -1), tr_abs.reshape(len(tr_abs), -1))
                t0 = cho.decide(te_abs.reshape(len(te_abs), -1))
                for p in lesion_p:
                    if p["file"] not in test:
                        continue
                    for c in contrasts:
                        t1 = float(cho.decide((p["absent"] + c * p["ell"]).reshape(1, -1))[0])
                        z[(p["file"], p["slice"], p["site_row"], p["site_col"], vol, c)] = per_lesion_z(t1, t0)
                for c in contrasts:
                    t1s = cho.decide(np.stack([p["absent"] + c * p["ell"] for p in lesion_p if p["file"] in test]).reshape(-1, roi * roi))
                    if len(t1s) > 1:
                        dprime_rows.append({"sequence": seq, "volume_mm3": vol, "contrast": c, "fold": k,
                                            "n_present": len(t1s), "n_absent": len(t0),
                                            "dprime_ref": float((t1s.mean() - t0.mean()) / np.sqrt(0.5 * (t1s.var(ddof=1) + t0.var(ddof=1))))})
    return z, pd.DataFrame(dprime_rows)


def save_patches(path, patches, roi):
    np.savez_compressed(path,
                        meta=json.dumps([{k: v for k, v in p.items() if k not in ("absent", "ell")} for p in patches]),
                        absent=np.stack([p["absent"] for p in patches]) if patches else np.zeros((0, roi, roi), np.float32),
                        ell=np.stack([p["ell"] if p["ell"] is not None else np.zeros((roi, roi), np.float32) for p in patches])
                        if patches else np.zeros((0, roi, roi), np.float32),
                        has_ell=np.array([p["ell"] is not None for p in patches]))


def load_patches(path):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    out = []
    for i, m in enumerate(meta):
        m = dict(m); m["absent"] = z["absent"][i]; m["ell"] = z["ell"][i] if bool(z["has_ell"][i]) else None
        out.append(m)
    return out


def apply_observer(rows, patches, roi, dprime_out):
    z, dprime = fit_reference_observer(patches, rows, roi)
    for r in rows:
        r["z_ref"] = z.get((r["file"], r["slice"], r["site_row"], r["site_col"], r["volume_mm3"], r["contrast"]), np.nan)
    dprime.to_csv(dprime_out, index=False)
    return rows, dprime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="outputs/bank_manifest.csv")
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="outputs/bank_characterization.csv")
    ap.add_argument("--patches-out", default="outputs/reference_patches.npz")
    ap.add_argument("--dprime-out", default="outputs/reference_dprime.csv")
    ap.add_argument("--maps-cache", default="/Volumes/T9/fastMRI_brain/cache/maps")
    ap.add_argument("--accelerations", default="4,6,8")
    ap.add_argument("--mask-types", default="equispaced")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cg-tol", type=float, default=1e-5)
    ap.add_argument("--cg-maxiter", type=int, default=80)
    ap.add_argument("--roi", type=int, default=32)
    ap.add_argument("--extra-absent", type=int, default=6)
    ap.add_argument("--limit-units", type=int, default=0, help="process at most N (file, slice) units (0 = all)")
    ap.add_argument("--observer-only", action="store_true",
                    help="skip pass 1; re-fit the reference observer from --out and --patches-out and refill z_ref")
    args = ap.parse_args()

    if args.observer_only:
        rows = pd.read_csv(args.out).to_dict("records")
        patches = load_patches(args.patches_out)
        rows, dprime = apply_observer(rows, patches, args.roi, args.dprime_out)
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"observer refit: z_ref available for {pd.DataFrame(rows)['z_ref'].notna().mean()*100:.0f}% of {len(rows)} rows; "
              f"{len(dprime)} d' rows -> {args.dprime_out}")
        return

    man = pd.read_csv(args.manifest)
    files = {os.path.splitext(f)[0]: os.path.join(dp, f) for dp, _, fs in os.walk(args.data) for f in fs if f.endswith(".h5")}
    done = set()
    if os.path.exists(args.out):
        prev = pd.read_csv(args.out)
        done = set(zip(prev["file"], prev["slice"]))
    jobs = []
    a = {"maps_cache": args.maps_cache, "accelerations": [int(x) for x in args.accelerations.split(",")],
         "mask_types": args.mask_types.split(","), "cg_tol": args.cg_tol, "cg_maxiter": args.cg_maxiter,
         "roi": args.roi, "extra_absent": args.extra_absent}
    for (stem, sl), g in man.groupby(["file", "slice"]):
        if stem not in files or (stem, sl) in done:
            continue
        jobs.append((files[stem], int(sl), g.to_dict("records"), a))
    if args.limit_units:
        jobs = jobs[:args.limit_units]
    print(f"{len(jobs)} (file, slice) units to process with {args.workers} workers; {len(done)} already done", flush=True)

    all_rows, all_patches = [], []
    t0 = time.time()
    with Pool(args.workers) as pool:
        for i, (rows, patches) in enumerate(pool.imap_unordered(process_unit, jobs), 1):
            all_rows.extend(rows); all_patches.extend(patches)
            if rows:
                pd.DataFrame(all_rows).to_csv(args.out + ".partial", index=False)
            print(f"[{time.time()-t0:6.0f}s] unit {i}/{len(jobs)}: {rows[0]['file'] if rows else '-'} slice {rows[0]['slice'] if rows else '-'} "
                  f"({len(rows)} rows, {rows[0]['unit_seconds'] if rows else 0}s)", flush=True)

    # merge with any earlier rows/patches (resumed run), then persist BEFORE the observer fit
    if os.path.exists(args.out):
        all_rows = pd.read_csv(args.out).to_dict("records") + all_rows
    if os.path.exists(args.patches_out):
        all_patches = load_patches(args.patches_out) + all_patches
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    pd.DataFrame(all_rows).to_csv(args.out, index=False)
    save_patches(args.patches_out, all_patches, args.roi)
    if os.path.exists(args.out + ".partial"):
        os.remove(args.out + ".partial")
    print(f"pass 1 saved: {len(all_rows)} rows -> {args.out}; {len(all_patches)} patches -> {args.patches_out}", flush=True)

    # pass 2: reference observer across files (re-runnable with --observer-only)
    all_rows, dprime = apply_observer(all_rows, all_patches, args.roi, args.dprime_out)
    df = pd.DataFrame(all_rows)
    df.to_csv(args.out, index=False)
    print(f"wrote {len(df)} rows -> {args.out}; {len(dprime)} d' rows -> {args.dprime_out}; "
          f"z_ref available for {df['z_ref'].notna().mean()*100:.0f}% of rows")


if __name__ == "__main__":
    main()
