#!/usr/bin/env python3
"""Reference detectability z_ref over an extended volume/contrast grid (CPU, cheap).

Uses the cached ESPIRiT maps and the manifest sites; inserts lesions into the
fully sampled reference only (no undersampling, no reconstruction, no CG),
extracts patches, and cross-fits the CHO by file exactly as 02b does. Purpose:
inform the pre-registration addendum on which lesion volumes and contrasts
are detectable in the reference on this data (partial volume: slices are
5-7.5 mm thick). No threshold is applied; the fraction above z_det is
reported for information only.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.lesion.insert import lesion_perturbation, slice_fill_fraction, sphere_radius_mm  # noqa: E402
from csa.observer.cho import CHO, extract_roi, laguerre_gauss_channels, per_lesion_z  # noqa: E402
from csa.physics.operator import ifft2c, sense_combine  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="outputs/annotated_val0/bank_manifest.csv")
    ap.add_argument("--data", default="/Volumes/T9/fastMRI_brain/extracted/multicoil_val")
    ap.add_argument("--maps-cache", default="/Volumes/T9/fastMRI_brain/cache/maps")
    ap.add_argument("--volumes", default="10,27,64,100,200,500")
    ap.add_argument("--contrasts-pos", default="0.24,0.5,1.0,1.5")
    ap.add_argument("--contrasts-neg", default="-0.3,-0.6,-1.0")
    ap.add_argument("--roi", type=int, default=32)
    ap.add_argument("--extra-absent", type=int, default=6)
    ap.add_argument("--thresholds", default="configs/thresholds.yaml")
    ap.add_argument("--out", default="outputs/annotated_val0/explore_reference_detectability.csv")
    ap.add_argument("--no-fill", action="store_true", help="also report without the partial-volume fill factor (thin-slice equivalent)")
    args = ap.parse_args()
    z_det = yaml.safe_load(open(args.thresholds)).get("z_det") or 3.0
    man = pd.read_csv(args.manifest)
    files = {os.path.splitext(f)[0]: os.path.join(dp, f) for dp, _, fs in os.walk(args.data) for f in fs if f.endswith(".h5")}
    vols = [float(v) for v in args.volumes.split(",")]
    roi = args.roi
    fills = [True, False] if args.no_fill else [True]
    patches = []  # dicts with absent patch, ell patch (unit contrast), meta
    t0 = time.time()
    rng = np.random.default_rng(0)
    for (stem, sl), g in man.groupby(["file", "slice"]):
        if stem not in files:
            continue
        with h5py.File(files[stem], "r") as f:
            ks = np.asarray(f["kspace"][sl]).astype(np.complex128)
        cp = os.path.join(args.maps_cache, f"{stem}_s{sl:02d}.npz")
        if not os.path.exists(cp):
            continue
        maps = np.load(cp)["maps"].astype(np.complex128)
        X = sense_combine(ifft2c(ks), maps); ssq = (np.abs(maps) ** 2).sum(0); support = ssq > 0.5
        sx = float(g["spacing_x_mm"].iloc[0] or 0.7); thick = float(g["thickness_mm"].iloc[0] or 5.0)
        seq = g["sequence"].iloc[0]
        sites = g.drop_duplicates(["site_row", "site_col"])[["site_row", "site_col"]].values
        for (row, col) in sites:
            site = (int(row), int(col))
            absent = extract_roi(X, site, roi).astype(np.float32)
            for v in vols:
                rmm = sphere_radius_mm(v); rpx = rmm / sx
                for use_fill in fills:
                    fill = slice_fill_fraction(rmm, thick) if use_fill else 1.0
                    try:
                        ell_u = lesion_perturbation(X, site, radius_px=rpx, contrast=1.0, fill=fill)
                    except ValueError:
                        continue
                    patches.append({"file": stem, "slice": sl, "sequence": seq, "site": site, "volume_mm3": v, "radius_px": rpx,
                                    "fill": fill, "with_fill": use_fill, "absent": absent, "ell": extract_roi(ell_u, site, roi).astype(np.float32)})
        a = np.abs(X); cand = np.argwhere(support & (a > 0.25 * a[support].max()))
        for k in range(args.extra_absent):
            r_, c_ = cand[rng.integers(len(cand))]
            patches.append({"file": stem, "slice": sl, "sequence": seq, "site": (int(r_), int(c_)), "volume_mm3": -1.0, "radius_px": 0.0,
                            "fill": 1.0, "with_fill": True, "absent": extract_roi(X, (r_, c_), roi).astype(np.float32), "ell": None})
    print(f"patches: {len(patches)} from {man.groupby(['file','slice']).ngroups} slices in {time.time()-t0:.0f}s", flush=True)

    rows = []
    for seq in sorted({p["sequence"] for p in patches}):
        files_seq = sorted({p["file"] for p in patches if p["sequence"] == seq})
        if len(files_seq) < 2:
            continue
        rng.shuffle(files_seq); folds = [set(files_seq[::2]), set(files_seq[1::2])]
        absent_all = [p for p in patches if p["sequence"] == seq]
        contrasts = [float(c) for c in (args.contrasts_pos if seq in ("AXFLAIR", "AXT2") else args.contrasts_neg).split(",")]
        for v in vols:
            for use_fill in fills:
                lesion_p = [p for p in patches if p["sequence"] == seq and p["volume_mm3"] == v and p["with_fill"] == use_fill]
                if not lesion_p:
                    continue
                rpx = float(np.median([p["radius_px"] for p in lesion_p]))
                U = laguerre_gauss_channels(roi, 6, a_px=float(np.clip(3.0 * rpx, 4.0, roi / 2)))
                c_fit = contrasts[len(contrasts) // 2]
                for k in (0, 1):
                    train, test = folds[k], folds[1 - k]
                    tr_abs = np.stack([p["absent"] for p in absent_all if p["file"] in train])
                    tr_pres = np.stack([p["absent"] + c_fit * p["ell"] for p in lesion_p if p["file"] in train])
                    te_abs = np.stack([p["absent"] for p in absent_all if p["file"] in test])
                    if len(tr_abs) < 8 or len(tr_pres) < 4 or len(te_abs) < 8:
                        continue
                    cho = CHO(U).fit(tr_pres.reshape(len(tr_pres), -1), tr_abs.reshape(len(tr_abs), -1))
                    t0s = cho.decide(te_abs.reshape(len(te_abs), -1))
                    for c in contrasts:
                        te_pres = [p for p in lesion_p if p["file"] in test]
                        t1s = cho.decide(np.stack([p["absent"] + c * p["ell"] for p in te_pres]).reshape(len(te_pres), -1))
                        zs = np.array([per_lesion_z(t, t0s) for t in t1s])
                        dprime = float((t1s.mean() - t0s.mean()) / np.sqrt(0.5 * (t1s.var(ddof=1) + t0s.var(ddof=1))))
                        rows.append({"sequence": seq, "volume_mm3": v, "with_fill": use_fill, "fill": float(np.median([p["fill"] for p in lesion_p])),
                                     "contrast": c, "fold": k, "n": len(zs), "z_median": float(np.median(zs)), "z_q25": float(np.quantile(zs, .25)),
                                     "z_q75": float(np.quantile(zs, .75)), "frac_ge_zdet": float((zs >= z_det).mean()), "dprime": dprime})
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    summ = df.groupby(["sequence", "with_fill", "volume_mm3", "contrast"])[["fill", "z_median", "frac_ge_zdet", "dprime"]].mean().round(2)
    print(summ.to_string())
    print(f"\nz_det = {z_det}; 'with_fill=False' rows are the thin-slice equivalent (no partial-volume dilution)")


if __name__ == "__main__":
    main()
