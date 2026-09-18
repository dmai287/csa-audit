#!/usr/bin/env python3
"""Build a small, self-contained pack of counterfactual test cases.

Each case is one real fastMRI slice with one lesion inserted through the real
forward operator, supplied in two forms:

  preview PNGs      the fully sampled reference, the same with the lesion, a
                    zero-filled reconstruction of each arm, and the lesion
                    transfer, for looking at;
  case .npz         everything another reconstruction model needs to run the
                    same comparison: the sampled k-space lines of both arms,
                    the sampling mask, the ESPIRiT coil maps, the reference
                    images and the lesion, plus the audit's own statistics.

Only the sampled lines are stored, not the zero-filled full grid, so a case is
a few megabytes rather than tens; `load_case()` in the pack's README rebuilds
the full array. Data are fastMRI and are covered by its data-use agreement.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from csa.lesion.insert import lesion_perturbation  # noqa: E402
from csa.metrics.image import center_crop  # noqa: E402
from csa.physics.masks import default_center_fraction, equispaced_mask  # noqa: E402
from csa.physics.operator import SenseOperator, ifft2c, sense_combine  # noqa: E402


def save_png(path, img, title, vmax=None, cmap="gray"):
    a = np.abs(img)
    vmax = vmax if vmax is not None else np.percentile(a, 99.5)
    fig, ax = plt.subplots(figsize=(4.2, 4.2 * a.shape[0] / max(a.shape[1], 1)))
    ax.imshow(a, cmap=cmap, vmin=0, vmax=vmax); ax.set_title(title, fontsize=8); ax.axis("off")
    fig.tight_layout(pad=0.2); fig.savefig(path, dpi=150); plt.close(fig)
    return vmax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", default="outputs/sample_pack_selection.csv")
    ap.add_argument("--data", default="/Volumes/T9/fastMRI_brain/extracted/multicoil_val")
    ap.add_argument("--maps-cache", default="/Volumes/T9/fastMRI_brain/cache/maps")
    ap.add_argument("--out", default="outputs/sample_pack")
    ap.add_argument("--acceleration", type=int, default=8)
    ap.add_argument("--crop", type=int, default=320)
    ap.add_argument("--calib", type=int, default=32, help="fully sampled centre lines shipped so maps can be recomputed")
    args = ap.parse_args()

    sel = pd.read_csv(args.selection).rename(columns={"lambda": "lam"})  # "lambda" is a keyword; itertuples renames it
    os.makedirs(args.out, exist_ok=True)
    index = []
    for i, r in enumerate(sel.itertuples(), 1):
        stem, sl = r.file, int(r.slice)
        cid = f"case{i:02d}_{r.sequence}_{int(r.volume_mm3)}mm3_c{r.contrast:+.2f}".replace("+", "p").replace("-", "m").replace(".", "")
        cdir = os.path.join(args.out, cid); os.makedirs(cdir, exist_ok=True)
        with h5py.File(os.path.join(args.data, stem + ".h5"), "r") as f:
            ks = np.asarray(f["kspace"][sl]).astype(np.complex128)
        maps = np.load(os.path.join(args.maps_cache, f"{stem}_s{sl:02d}.npz"))["maps"].astype(np.complex128)
        X = sense_combine(ifft2c(ks), maps)
        mask = equispaced_mask(ks.shape[-1], args.acceleration, default_center_fraction(args.acceleration))
        op = SenseOperator(maps, mask)
        ell = lesion_perturbation(X, (int(r.site_row), int(r.site_col)), radius_px=float(r.radius_px),
                                  contrast=float(r.contrast), fill=float(r.fill))
        X_cf = X + ell
        y, y_cf = op.undersample(ks), op.undersample(ks) + op.forward(ell)
        calib_lo = (ks.shape[-1] - args.calib) // 2
        zf, zf_cf = op.adjoint(y), op.adjoint(y_cf)

        c = lambda a: center_crop(a, args.crop)
        vmax = save_png(os.path.join(cdir, "1_reference.png"), c(X), "fully sampled reference (lesion absent)")
        save_png(os.path.join(cdir, "2_reference_with_lesion.png"), c(X_cf), f"reference + lesion ({r.volume_mm3:g} mm3, contrast {r.contrast:+.2f})", vmax)
        save_png(os.path.join(cdir, "3_zerofilled_lesion_absent.png"), c(zf), f"zero-filled, R={args.acceleration}, lesion absent", vmax)
        save_png(os.path.join(cdir, "4_zerofilled_lesion_present.png"), c(zf_cf), f"zero-filled, R={args.acceleration}, lesion present", vmax)
        save_png(os.path.join(cdir, "5_lesion_transfer.png"), c(zf_cf - zf), "lesion transfer (present minus absent)", None, "magma")
        save_png(os.path.join(cdir, "6_lesion_truth.png"), c(ell), "the inserted lesion itself", None, "magma")

        np.savez_compressed(
            os.path.join(cdir, "case.npz"),
            kspace_sampled_lesion_absent=y[:, :, mask].astype(np.complex64),
            kspace_sampled_lesion_present=y_cf[:, :, mask].astype(np.complex64),
            mask=mask.astype(bool),
            calibration_kspace=ks[:, :, calib_lo:calib_lo + args.calib].astype(np.complex64),
            calibration_first_line=np.int32(calib_lo),
            reference_lesion_absent=X.astype(np.complex64),
            reference_lesion_present=X_cf.astype(np.complex64),
            lesion=ell.astype(np.complex64),
            meta=json.dumps({
                "case_id": cid, "fastmri_file": stem, "slice": sl, "sequence": r.sequence,
                "n_coils": int(ks.shape[0]), "image_shape": list(ks.shape[-2:]),
                "acceleration": args.acceleration, "mask_type": "equispaced",
                "sampled_lines": int(mask.sum()), "total_lines": int(mask.size),
                "calibration_lines": int(args.calib),
                "lesion": {"volume_mm3": float(r.volume_mm3), "contrast": float(r.contrast),
                           "radius_px": float(r.radius_px), "fill_fraction": float(r.fill),
                           "site_row": int(r.site_row), "site_col": int(r.site_col)},
                "audit_statistics": {"measurement_gain_kappa2": float(r.kappa2),
                                     "measured_fraction_mu_lambda": float(r.mu_lambda),
                                     "reference_detectability_z_ref": float(r.z_ref),
                                     "lambda": float(r.lam)},
                "note": "kspace_sampled_* holds only the sampled phase-encode lines; scatter them "
                        "back with the mask to rebuild the full grid. Coil maps are not shipped: "
                        "recompute them from calibration_kspace with ESPIRiT, as the audit does, so "
                        "your maps are produced the same way ours were."}))
        size = sum(os.path.getsize(os.path.join(cdir, f)) for f in os.listdir(cdir)) / 1e6
        index.append({"case_id": cid, "fastmri_file": stem, "slice": sl, "sequence": r.sequence,
                      "volume_mm3": r.volume_mm3, "contrast": r.contrast, "n_coils": int(ks.shape[0]),
                      "mu_lambda": round(float(r.mu_lambda), 3), "kappa2": round(float(r.kappa2), 3),
                      "z_ref": round(float(r.z_ref), 2), "size_MB": round(size, 1)})
        print(f"[{i}/{len(sel)}] {cid}: {size:.1f} MB", flush=True)
    pd.DataFrame(index).to_csv(os.path.join(args.out, "index.csv"), index=False)
    print("\n" + pd.DataFrame(index).to_string(index=False))
    print(f"\ntotal {sum(r['size_MB'] for r in index):.1f} MB in {args.out}")


if __name__ == "__main__":
    main()
