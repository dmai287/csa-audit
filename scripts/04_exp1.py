#!/usr/bin/env python3
"""Experiment 1: counterfactual lesion stress test, per pair, with pluggable reconstructors.

Runs Algorithm 1 for every planned insertion in the bank manifest and returns
one row per (pair, model, acceleration): both arms reconstructed under a
common seed, the measured and null-space transfer of the lesion, image scores,
ROI measures, and the observer readout on the reconstruction with two-fold
cross-fitting by file. Erasure and silent-erasure flags are computed from the
pre-registered thresholds in configs/thresholds.yaml.

Which reconstructors run is a command-line choice. With the classical
reconstructors (`zero_filled`, `cg_sense`) this is a CPU-only DRY RUN of the
full machinery and a "no prior" reference; the confirmatory experiment adds
the learned models on a GPU and is subject to the pre-registration's order
(OSF deposit, then pilot, then this). Output files carry the run label in
their name so the two are never confused.

Linear reconstructors: the counterfactual arm is reconstructed once at unit
contrast and scaled, which is exact for a linear map and cuts the cost by the
number of contrast levels.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from multiprocessing import Pool

import h5py
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.audit.erasure import is_erased, within_central_interval  # noqa: E402
from csa.lesion.insert import lesion_perturbation, soft_disc  # noqa: E402
from csa.metrics.image import annulus, center_crop, psnr, roi_cnr, roi_ssim, ssim  # noqa: E402
from csa.nullspace.projectors import measurement_gain, null_filter, real_inner  # noqa: E402
from csa.observer.cho import CHO, extract_roi, laguerre_gauss_channels, per_lesion_z  # noqa: E402
from csa.physics.masks import default_center_fraction, equispaced_mask, random_mask  # noqa: E402
from csa.physics.operator import SenseOperator, ifft2c, sense_combine  # noqa: E402
from csa.recon.base import CGSense, ZeroFilled  # noqa: E402

LINEAR = {"zero_filled", "cg_sense"}


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def make_reconstructor(name, lam, cg_tol, cg_maxiter):
    if name == "zero_filled":
        return ZeroFilled()
    if name == "cg_sense":
        return CGSense(lam=lam, tol=cg_tol, maxiter=cg_maxiter)
    raise ValueError(f"reconstructor '{name}' needs a GPU adapter; see docs/MODELS.md")


def process_unit(job):
    (path, sl, model, R, mt, rows_spec, char_rows, args) = job
    t_unit = time.time()
    with h5py.File(path, "r") as f:
        ks = np.asarray(f["kspace"][sl]).astype(np.complex128)
    stem = os.path.splitext(os.path.basename(path))[0]
    maps = np.load(os.path.join(args["maps_cache"], f"{stem}_s{sl:02d}.npz"))["maps"].astype(np.complex128)
    X = sense_combine(ifft2c(ks), maps)
    ssq = (np.abs(maps) ** 2).sum(0); support = ssq > 0.5
    lam = float(np.median([c["lambda"] for c in char_rows]))
    cf = default_center_fraction(R)
    mask = equispaced_mask(ks.shape[-1], R, cf) if mt == "equispaced" else random_mask(ks.shape[-1], R, cf, seed=0)
    op = SenseOperator(maps, mask)
    rec = make_reconstructor(model, lam, args["cg_tol"], args["cg_maxiter"])
    seed = args["seed"]
    y = op.undersample(ks)
    Xh = rec(y, mask, maps, seed)
    res_fact = op.relative_residual(Xh, y)
    crop = (lambda a: center_crop(a, args["metric_crop"])) if args.get("metric_crop") else (lambda a: a)
    psnr_fact, ssim_fact = psnr(crop(X), crop(Xh)), ssim(crop(X), crop(Xh))
    roi = args["roi"]
    rng = np.random.default_rng((hash((stem, sl, model, R)) % (2 ** 31)))
    char = {(int(c["site_row"]), int(c["site_col"]), float(c["volume_mm3"])): c for c in char_rows}

    seen = {}
    for spec in rows_spec:
        key = (int(spec["site_row"]), int(spec["site_col"]), float(spec["volume_mm3"]))
        seen.setdefault(key, []).append(float(spec["contrast"]))
    out_rows, patches = [], []
    n_zero_gain = 0
    for (row, col, vol), contrasts in seen.items():
        c_info = char.get((row, col, vol))
        if c_info is None:
            continue
        rpx, fill = float(c_info["radius_px"]), float(c_info["fill"])
        site = (row, col)
        try:
            ell_u = lesion_perturbation(X, site, radius_px=rpx, contrast=1.0, fill=fill)
        except ValueError:
            continue
        if not np.any(ell_u):
            continue
        Al_u = op.forward(ell_u)
        if real_inner(Al_u, Al_u) == 0.0:
            n_zero_gain += 1  # site outside the coil-map support: nothing reaches k-space
            continue
        disc = soft_disc(X.shape, site, rpx, taper_px=0.0) > 0.5
        ring = annulus(X.shape, site, rpx + 2, rpx + 8) & support
        Ql_u = null_filter(op, ell_u, lam, tol=args["cg_tol"], maxiter=args["cg_maxiter"])
        absent_patch = extract_roi(Xh, site, roi).astype(np.float32)
        linear = model in LINEAR
        if linear:
            dX_u = rec(y + op.forward(ell_u), mask, maps, seed) - Xh
            QdX_u = null_filter(op, dX_u, lam, tol=args["cg_tol"], maxiter=args["cg_maxiter"])
        for c in contrasts:
            if c == 0.0:
                continue  # a zero-contrast lesion has no transfer to measure
            t0 = time.time()
            ell = c * ell_u
            if linear:
                dX = c * dX_u; QdX = c * QdX_u
            else:
                Xh_cf_full = rec(y + op.forward(ell), mask, maps, seed)
                dX = Xh_cf_full - Xh
                QdX = null_filter(op, dX, lam, tol=args["cg_tol"], maxiter=args["cg_maxiter"])
            Xh_cf = Xh + dX
            X_cf = X + ell
            y_cf = y + op.forward(ell)
            Al, Ad = c * Al_u, op.forward(dX)
            t_R = real_inner(Ad, Al) / real_inner(Al, Al)
            Ql = c * Ql_u
            t_N = real_inner(QdX, Ql) / real_inner(Ql, Ql)
            present_patch = extract_roi(Xh_cf, site, roi).astype(np.float32)
            patches.append({"file": stem, "slice": sl, "model": model, "acceleration": R, "mask_type": mt,
                            "sequence": rows_spec[0]["sequence"], "site_row": row, "site_col": col, "volume_mm3": vol,
                            "contrast": c, "radius_px": rpx, "kind": "present", "patch": present_patch})
            out_rows.append({
                "file": stem, "slice": sl, "sequence": rows_spec[0]["sequence"], "site_row": row, "site_col": col,
                "volume_mm3": vol, "contrast": c, "radius_px": rpx, "fill": fill, "model": model,
                "acceleration": R, "mask_type": mt, "seed": seed, "lambda": lam,
                "kappa": measurement_gain(op, ell), "mu_lambda": float(c_info.get("mu_lambda", np.nan)),
                "z_ref": float(c_info.get("z_ref", np.nan)) if "z_ref" in c_info else np.nan,
                "residual_factual": res_fact, "residual_cf": op.relative_residual(Xh_cf, y_cf),
                "t_R": float(t_R), "t_N": float(t_N),
                "psnr_factual": psnr_fact, "psnr_cf": psnr(crop(X_cf), crop(Xh_cf)), "ssim_factual": ssim_fact, "ssim_cf": ssim(crop(X_cf), crop(Xh_cf)),
                "metric_crop": args.get("metric_crop") or 0,
                "roi_cnr_cf": roi_cnr(Xh_cf, disc, ring) if ring.any() else np.nan,
                "roi_cnr_ref": roi_cnr(X_cf, disc, ring) if ring.any() else np.nan,
                "roi_ssim_cf": roi_ssim(X_cf, Xh_cf, site, roi),
                "lesion_energy": float(real_inner(ell, ell)), "pair_seconds": round(time.time() - t0, 1)})
        patches.append({"file": stem, "slice": sl, "model": model, "acceleration": R, "mask_type": mt,
                        "sequence": rows_spec[0]["sequence"], "site_row": row, "site_col": col, "volume_mm3": vol,
                        "contrast": 0.0, "radius_px": rpx, "kind": "absent", "patch": absent_patch})
    a = np.abs(Xh); cand = np.argwhere(support & (a > 0.25 * a[support].max()))
    for k in range(args["extra_absent"]):
        r_, c_ = cand[rng.integers(len(cand))]
        patches.append({"file": stem, "slice": sl, "model": model, "acceleration": R, "mask_type": mt,
                        "sequence": rows_spec[0]["sequence"], "site_row": int(r_), "site_col": int(c_), "volume_mm3": -1.0,
                        "contrast": 0.0, "radius_px": 0.0, "kind": "absent", "patch": extract_roi(Xh, (r_, c_), roi).astype(np.float32)})
    for rr in out_rows:
        rr["unit_seconds"] = round(time.time() - t_unit, 1); rr["sites_skipped_zero_gain"] = n_zero_gain
    return out_rows, patches, f"{stem}_s{sl:02d}_{model}_R{R}_{mt}"


def unit_path(unit_dir, key):
    return os.path.join(unit_dir, key + ".npz")


def save_unit(unit_dir, key, rows, patches, roi):
    os.makedirs(unit_dir, exist_ok=True)
    np.savez_compressed(unit_path(unit_dir, key), rows=json.dumps(rows),
                        meta=json.dumps([{k: v for k, v in p.items() if k != "patch"} for p in patches]),
                        patch=np.stack([p["patch"] for p in patches]) if patches else np.zeros((0, roi, roi), np.float32))


def load_units(unit_dir):
    rows, patches = [], []
    for f in sorted(os.listdir(unit_dir)) if os.path.isdir(unit_dir) else []:
        if not f.endswith(".npz"):
            continue
        z = np.load(os.path.join(unit_dir, f), allow_pickle=False)
        rows.extend(json.loads(str(z["rows"])))
        for i, m in enumerate(json.loads(str(z["meta"]))):
            m = dict(m); m["patch"] = z["patch"][i]; patches.append(m)
    return rows, patches


def fit_recon_observers(patches, rows, roi, seed=0):
    """Two-fold cross-fitted CHO per (model, R, mask, sequence, volume) on reconstructed patches."""
    z = {}; dp_rows = []
    rng = np.random.default_rng(seed)
    meta = pd.DataFrame([{k: v for k, v in p.items() if k != "patch"} for p in patches])
    if meta.empty:
        return z, pd.DataFrame(dp_rows)
    for (model, R, mt, seq), gmeta in meta.groupby(["model", "acceleration", "mask_type", "sequence"]):
        files = sorted(gmeta["file"].unique())
        if len(files) < 2:
            continue
        rng.shuffle(files); folds = [set(files[::2]), set(files[1::2])]
        grp = [p for p in patches if (p["model"], p["acceleration"], p["mask_type"], p["sequence"]) == (model, R, mt, seq)]
        absent = [p for p in grp if p["kind"] == "absent"]
        for vol in sorted(v for v in gmeta["volume_mm3"].unique() if v > 0):
            pres = [p for p in grp if p["kind"] == "present" and p["volume_mm3"] == vol]
            if not pres:
                continue
            rpx = float(np.median([p["radius_px"] for p in pres]))
            U = laguerre_gauss_channels(roi, 6, a_px=float(np.clip(3.0 * rpx, 4.0, roi / 2)))
            contrasts = sorted({p["contrast"] for p in pres}); c_fit = contrasts[len(contrasts) // 2]
            for k in (0, 1):
                train, test = folds[k], folds[1 - k]
                tr_abs = [p["patch"] for p in absent if p["file"] in train]
                tr_pres = [p["patch"] for p in pres if p["file"] in train and p["contrast"] == c_fit]
                te_abs = [p["patch"] for p in absent if p["file"] in test]
                if len(tr_abs) < 8 or len(tr_pres) < 4 or len(te_abs) < 8:
                    continue
                cho = CHO(U).fit(np.stack(tr_pres).reshape(len(tr_pres), -1), np.stack(tr_abs).reshape(len(tr_abs), -1))
                t0 = cho.decide(np.stack(te_abs).reshape(len(te_abs), -1))
                for p in pres:
                    if p["file"] in test:
                        t1 = float(cho.decide(p["patch"].reshape(1, -1))[0])
                        z[(p["file"], p["slice"], model, R, mt, p["site_row"], p["site_col"], vol, p["contrast"])] = per_lesion_z(t1, t0)
                for c in contrasts:
                    t1s = cho.decide(np.stack([p["patch"] for p in pres if p["file"] in test and p["contrast"] == c]).reshape(-1, roi * roi)) \
                        if any(p["file"] in test and p["contrast"] == c for p in pres) else np.array([])
                    if len(t1s) > 1:
                        dp_rows.append({"model": model, "acceleration": R, "mask_type": mt, "sequence": seq, "volume_mm3": vol,
                                        "contrast": c, "fold": k, "n_present": len(t1s), "n_absent": len(t0),
                                        "dprime_recon": float((t1s.mean() - t0.mean()) / np.sqrt(0.5 * (t1s.var(ddof=1) + t0.var(ddof=1))))})
    return z, pd.DataFrame(dp_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--characterization", required=True, help="output of 02b (lambda, radius_px, fill, mu_lambda, z_ref)")
    ap.add_argument("--data", required=True)
    ap.add_argument("--maps-cache", default="/Volumes/T9/fastMRI_brain/cache/maps")
    ap.add_argument("--models", default="zero_filled,cg_sense")
    ap.add_argument("--accelerations", default="4,8")
    ap.add_argument("--mask-types", default="equispaced")
    ap.add_argument("--slices-per-file", type=int, default=0, help="use only the central N manifest slices per file (0 = all)")
    ap.add_argument("--thresholds", default="configs/thresholds.yaml")
    ap.add_argument("--label", default="dryrun_classical")
    ap.add_argument("--out-dir", default="outputs/annotated_val0")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cg-tol", type=float, default=1e-5)
    ap.add_argument("--cg-maxiter", type=int, default=80)
    ap.add_argument("--roi", type=int, default=32)
    ap.add_argument("--extra-absent", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--metric-crop", type=int, default=320,
                    help="global PSNR/SSIM on the central N x N crop as fastMRI reports them (0 = full grid)")
    ap.add_argument("--limit-units", type=int, default=0)
    ap.add_argument("--observer-only", action="store_true", help="skip pass 1; refit observers and erasure flags from saved rows and patches")
    args = ap.parse_args()

    out_csv = os.path.join(args.out_dir, f"exp1_{args.label}.csv")
    unit_dir = os.path.join(args.out_dir, f"exp1_{args.label}_units")
    dp_out = os.path.join(args.out_dir, f"exp1_{args.label}_dprime.csv")
    thr = yaml.safe_load(open(args.thresholds))
    if args.observer_only:
        rows, patches = load_units(unit_dir)
        df, dp = pass2(rows, patches, args, thr)
        df.to_csv(out_csv, index=False); dp.to_csv(dp_out, index=False)
        report(df, dp, args, out_csv, dp_out); return
    man = pd.read_csv(args.manifest)
    char = pd.read_csv(args.characterization)
    files = {os.path.splitext(f)[0]: os.path.join(dp, f) for dp, _, fs in os.walk(args.data) for f in fs if f.endswith(".h5")}
    a = {"maps_cache": args.maps_cache, "cg_tol": args.cg_tol, "cg_maxiter": args.cg_maxiter, "roi": args.roi,
         "extra_absent": args.extra_absent, "seed": args.seed, "metric_crop": args.metric_crop}
    jobs, n_done = [], 0
    for stem, g_file in man.groupby("file"):
        if stem not in files:
            continue
        slices = sorted(g_file["slice"].unique())
        if args.slices_per_file:
            lo = max((len(slices) - args.slices_per_file) // 2, 0); slices = slices[lo:lo + args.slices_per_file]
        for sl in slices:
            rows_spec = g_file[g_file["slice"] == sl].to_dict("records")
            char_rows = char[(char["file"] == stem) & (char["slice"] == sl)]
            if char_rows.empty:
                continue
            for model in args.models.split(","):
                for R in (int(x) for x in args.accelerations.split(",")):
                    for mt in args.mask_types.split(","):
                        if os.path.exists(unit_path(unit_dir, f"{stem}_s{int(sl):02d}_{model}_R{R}_{mt}")):
                            n_done += 1; continue
                        cr = char_rows[(char_rows["acceleration"] == R) & (char_rows["mask_type"] == mt)]
                        if cr.empty:
                            continue
                        jobs.append((files[stem], int(sl), model, R, mt, rows_spec, cr.to_dict("records"), a))
    if args.limit_units:
        jobs = jobs[:args.limit_units]
    print(f"[{args.label}] {len(jobs)} units (file, slice, model, R, mask) with {args.workers} workers; {n_done} already done", flush=True)
    if any(j[2] not in LINEAR for j in jobs):
        print("NOTE: learned models requested; make sure the pre-registration's order (OSF deposit, pilot, addendum) has been honoured.", flush=True)

    t0 = time.time()
    if jobs:
        with Pool(args.workers) as pool:
            for i, (rows, patches, key) in enumerate(pool.imap_unordered(process_unit, jobs), 1):
                save_unit(unit_dir, key, rows, patches, args.roi)
                print(f"[{time.time()-t0:6.0f}s] unit {i}/{len(jobs)}: {key} ({len(rows)} pairs, "
                      f"{rows[0]['unit_seconds'] if rows else 0}s, {rows[0]['sites_skipped_zero_gain'] if rows else '-'} sites skipped)", flush=True)

    all_rows, all_patches = load_units(unit_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    df, dp = pass2(all_rows, all_patches, args, thr)
    df.to_csv(out_csv, index=False); dp.to_csv(dp_out, index=False)
    report(df, dp, args, out_csv, dp_out)


def pass2(all_rows, all_patches, args, thr):
    """Observers on reconstructions, then erasure flags from the pre-registered thresholds."""
    z, dp = fit_recon_observers(all_patches, all_rows, args.roi)
    df = pd.DataFrame(all_rows)
    df["z_recon"] = [z.get((r.file, r.slice, r.model, r.acceleration, r.mask_type, r.site_row, r.site_col, r.volume_mm3, r.contrast), np.nan)
                     for r in df.itertuples()]
    z_det, z_miss = thr.get("z_det"), thr.get("z_miss"); level = float(thr.get("psnr_interval_level", 0.95))
    if z_det is None or z_miss is None:
        print("thresholds are null; erasure flags not computed"); df["erased"] = np.nan; df["silently_erased"] = np.nan
    else:
        df["erased"] = [bool(is_erased(zr, zc, z_det, z_miss)) if np.isfinite(zr) and np.isfinite(zc) else np.nan
                        for zr, zc in zip(df["z_ref"], df["z_recon"])]
        df["indeterminate"] = [bool(np.isfinite(zr) and np.isfinite(zc) and zr >= z_det and z_miss <= zc < z_det) for zr, zc in zip(df["z_ref"], df["z_recon"])]
        # lesion-free PSNR distribution: factual-arm PSNR across slices of the same model/R/mask
        silent = []
        for r in df.itertuples():
            pool_ = df[(df.model == r.model) & (df.acceleration == r.acceleration) & (df.mask_type == r.mask_type)].drop_duplicates(["file", "slice"])["psnr_factual"].values
            silent.append(bool(r.erased) and within_central_interval(r.psnr_cf, pool_, level) if r.erased is True or r.erased == 1.0 else (np.nan if r.erased != r.erased else False))
        df["silently_erased"] = silent
    df["git_commit"] = git_commit(); df["thresholds_sha256"] = hashlib.sha256(open(args.thresholds, "rb").read()).hexdigest()[:12]
    df["run_label"] = args.label
    return df, dp


def report(df, dp, args, out_csv, dp_out):
    print(f"[{args.label}] wrote {len(df)} pair rows -> {out_csv}; z_recon for {df['z_recon'].notna().mean()*100:.0f}%; {len(dp)} d' rows -> {dp_out}")
    if "erased" in df and df["erased"].notna().any():
        summ = df.dropna(subset=["erased"]).groupby(["model", "acceleration"]).agg(pairs=("erased", "size"), erased=("erased", "mean"), silent=("silently_erased", "mean"), t_N=("t_N", "mean"), t_R=("t_R", "mean")).round(3)
        print(f"[{args.label}] DRY-RUN summary (not a confirmatory result):\n{summ}")


if __name__ == "__main__":
    main()
