#!/usr/bin/env python3
"""Build the synthetic lesion-bank manifest and per-sequence appearance statistics.

Two outputs:
  outputs/appearance_stats.csv   -- one row per fastMRI+ bounding box: sequence,
                                     bbox magnitude, local background magnitude,
                                     contrast ratio. Read once from real data so
                                     configs/bank.yaml's contrast levels are set
                                     from evidence, not guessed.
  outputs/bank_manifest.csv      -- one row per planned synthetic insertion:
                                     host file/slice, site, volume, contrast,
                                     drawn from lesion-free slices only.

Coverage note: fastMRI+ brain annotations exist only for AXFLAIR, AXT1 and
AXT1POST; there are no AXT2 or hemorrhage-labelled rows in this release. T2
hosts therefore get the synthetic (Case A) bank only, with no Case B
real-lesion counterpart, and this is recorded as a limitation.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import h5py
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.metrics.image import annulus  # noqa: E402
from csa.physics.operator import ifft2c  # noqa: E402
from csa.io.fastmri import pixel_spacing_mm  # noqa: E402

SEQ_RE = re.compile(r"file_brain_(AX[A-Z0-9]+)_")


def sequence_of(fname: str) -> str:
    m = SEQ_RE.match(fname)
    return m.group(1) if m else "UNKNOWN"


def rss_slice(kspace_slice: np.ndarray) -> np.ndarray:
    """Root-sum-of-squares magnitude image from one fully sampled multi-coil slice.

    Used only to measure lesion appearance statistics from real annotations;
    the audit itself uses ESPIRiT + sensitivity combination (Section 3.2), not RSS.
    """
    coil_images = ifft2c(kspace_slice)
    return np.sqrt((np.abs(coil_images) ** 2).sum(axis=0))


def load_annotations(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[(df["study_level"] == "No") & (df["width"] > 0) & (df["height"] > 0)].copy()
    df["sequence"] = df["file"].map(sequence_of)
    return df


def index_local_files(data_dir: str) -> dict:
    """Map fastMRI file stem -> local .h5 path, across whichever subdirectories exist."""
    idx = {}
    for path in glob.glob(os.path.join(data_dir, "**", "file_brain_*.h5"), recursive=True):
        idx[os.path.splitext(os.path.basename(path))[0]] = path
    return idx


FOCAL_LABELS = {"Nonspecific white matter lesion", "Nonspecific lesion", "Lacunar infarct", "Mass",
                "Extra-axial mass", "Encephalomalacia", "Small vessel chronic white matter ischemic change"}
MAX_FOCAL_SIDE_PX = 25   # about 17 mm at 0.7 mm/px; larger boxes are structural, not focal


def appearance_stats(annotations: pd.DataFrame, file_index: dict, ring_px: int = 8) -> pd.DataFrame:
    """For each annotated bounding box present locally: RSS magnitude statistics of
    the box against an annulus around it.

    Two contrast measures are recorded. `contrast_mean` is the box mean against
    the annulus mean, which a loose bounding box dilutes toward zero.
    `contrast_core` is the signed deviation of the pixel at the 90th percentile
    of |x - bg| / bg inside the box, i.e. the lesion core, and is the statistic
    the bank uses. `focal` marks boxes whose label is a focal lesion and whose
    longest side is at most MAX_FOCAL_SIDE_PX."""
    rows = []
    open_files = {}
    try:
        for _, r in annotations.iterrows():
            path = file_index.get(r["file"])
            if path is None:
                continue
            if path not in open_files:
                open_files[path] = h5py.File(path, "r")
            f = open_files[path]
            sl = int(r["slice"])
            if sl >= f["kspace"].shape[0]:
                continue
            img = rss_slice(np.asarray(f["kspace"][sl]))
            x, y, w, h = int(r["x"]), int(r["y"]), int(r["width"]), int(r["height"])
            cy, cx = y + h / 2.0, x + w / 2.0
            box = np.zeros(img.shape, dtype=bool)
            box[max(y, 0):y + h, max(x, 0):x + w] = True
            if not box.any():
                continue
            radius = max(w, h) / 2.0
            ring = annulus(img.shape, (cy, cx), radius + 2, radius + 2 + ring_px)
            if not ring.any():
                continue
            box_mag, bg_mag = float(img[box].mean()), float(img[ring].mean())
            if bg_mag <= 0:
                continue
            dev = (img[box] - bg_mag) / bg_mag
            k = int(np.argsort(np.abs(dev))[int(0.9 * (dev.size - 1))])
            rows.append({"file": r["file"], "slice": sl, "sequence": r["sequence"],
                        "label": r["label"], "width": w, "height": h,
                        "focal": bool(r["label"] in FOCAL_LABELS and max(w, h) <= MAX_FOCAL_SIDE_PX),
                        "box_mag": box_mag, "background_mag": bg_mag,
                        "contrast_mean": (box_mag - bg_mag) / bg_mag,
                        "contrast_core": float(dev[k])})
    finally:
        for f in open_files.values():
            f.close()
    cols = ["file", "slice", "sequence", "label", "width", "height", "focal", "box_mag", "background_mag",
            "contrast_mean", "contrast_core"]
    return pd.DataFrame(rows, columns=cols)


def annotated_slices(annotations: pd.DataFrame) -> set:
    """(file, slice) pairs to exclude from insertion, regardless of sequence."""
    return set(zip(annotations["file"], annotations["slice"].astype(int)))


def build_manifest(file_index: dict, excluded: set, bank_cfg: dict, contrast_by_seq: dict,
                   seed: int = 0, middle_slices: int = 0, maps_cache: str = None) -> pd.DataFrame:
    """Sample (file, slice, site, volume, contrast) from lesion-free slices.

    The three contrast levels share one site per (slice, volume), so the
    model-independent quantities (measurement gain and fraction) are computed
    once per site and reused; contrast only scales the perturbation.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for stem, path in sorted(file_index.items()):
        seq = sequence_of(stem)
        contrasts = contrast_by_seq.get(seq)
        if not contrasts:
            continue  # no evidence-based contrast for this sequence yet
        with h5py.File(path, "r") as f:
            n_slices = f["kspace"].shape[0]
            shape = f["kspace"].shape[-2:]
            header = f["ismrmrd_header"][()]
        header = header.decode() if isinstance(header, bytes) else str(header)
        try:
            spacing = pixel_spacing_mm(header)
        except Exception:
            spacing = {}
        slices = range(n_slices)
        if middle_slices > 0:
            lo = max((n_slices - middle_slices) // 2, 0)
            slices = range(lo, min(lo + middle_slices, n_slices))
        for sl in slices:
            if (stem, sl) in excluded:
                continue
            with h5py.File(path, "r") as f:
                img = rss_slice(np.asarray(f["kspace"][sl]))
            a = img / img.max()
            ok = (a > 0.25) & (a < 0.9)
            cp = os.path.join(maps_cache, f"{stem}_s{sl:02d}.npz") if maps_cache else None
            if cp and os.path.exists(cp):
                m = np.load(cp)["maps"]; ok &= (np.abs(m) ** 2).sum(0) > 0.5  # inside the ESPIRiT support
            candidates = np.argwhere(ok)
            if len(candidates) == 0:
                continue
            for volume in bank_cfg["volumes_mm3"]:
                row, col = (int(v) for v in candidates[rng.integers(len(candidates))])
                for contrast in contrasts:
                    rows.append({"file": stem, "slice": sl, "sequence": seq,
                                "site_row": row, "site_col": col,
                                "volume_mm3": volume, "contrast": contrast,
                                "spacing_x_mm": spacing.get("x"), "spacing_y_mm": spacing.get("y"),
                                "thickness_mm": spacing.get("z"), "n_rows": shape[0], "n_cols": shape[1]})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations", default="/Volumes/T9/fastMRI_brain/fastmri_plus/Annotations/brain.csv")
    ap.add_argument("--data", required=True, help="directory of extracted fastMRI brain .h5 files")
    ap.add_argument("--config", default="configs/bank.yaml")
    ap.add_argument("--stats-out", default="outputs/appearance_stats.csv")
    ap.add_argument("--manifest-out", default="outputs/bank_manifest.csv")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--middle-slices", type=int, default=0, help="restrict to the central N slices per volume (0 = all)")
    ap.add_argument("--maps-cache", default="/Volumes/T9/fastMRI_brain/cache/maps",
                    help="if a cached ESPIRiT map exists for a slice, insertion sites are restricted to its support")
    ap.add_argument("--contrast-levels", default="",
                    help="force levels per sequence, e.g. 'AXFLAIR=0.24,1.0,1.5;AXT1=-0.21,-0.6,-1.0' "
                         "(recorded as contrast_source=forced; use for a pre-registration amendment)")
    ap.add_argument("--default-contrasts", default="",
                    help="comma-separated contrast levels for sequences with no fastMRI+ statistics "
                         "(e.g. AXT2, which has no annotations); recorded as contrast_source=default")
    args = ap.parse_args()

    bank_cfg = yaml.safe_load(open(args.config))
    annotations = load_annotations(args.annotations)
    file_index = index_local_files(args.data)
    print(f"{len(file_index)} local files indexed; {len(annotations)} slice-level annotations total")

    stats = appearance_stats(annotations, file_index)
    os.makedirs(os.path.dirname(args.stats_out) or ".", exist_ok=True)
    stats.to_csv(args.stats_out, index=False)
    print(f"appearance stats measured on {len(stats)} locally-present annotated boxes -> {args.stats_out}")
    if stats.empty:
        print("WARNING: none of the annotated files are present locally; "
              "contrast levels cannot be derived yet. Extract more archives, or "
              "pass --data pointing at a directory that includes annotated files "
              "(AXFLAIR/AXT1/AXT1POST volumes -- see brain_file_list.csv).")

    contrast_by_seq, contrast_source = {}, {}
    signs = bank_cfg.get("contrast_signs", {})
    floor = float(bank_cfg.get("contrast_floor", 0.1))
    min_focal = int(bank_cfg.get("min_focal_boxes", 10))
    for seq, g in (stats.groupby("sequence") if len(stats) else []):
        f = g[g["focal"]]
        print(f"  {seq}: {len(g)} boxes, {len(f)} focal; core |contrast| quartiles over focal boxes: "
              + (", ".join(f"{v:.2f}" for v in f["contrast_core"].abs().quantile([0.25, 0.5, 0.75])) if len(f) else "n/a"))
        if len(f) < min_focal:
            continue
        q1, q2, q3 = f["contrast_core"].abs().quantile([0.25, 0.5, 0.75])
        sign = signs.get(seq, signs.get(seq.replace("AX", ""), [1]))
        sign = float(sign[0]) if isinstance(sign, (list, tuple)) else float(sign)
        levels = sorted({round(sign * max(q1, floor), 2), round(sign * max(q2, floor), 2), round(sign * max(q3, floor), 2)})
        source = "fastmri_plus_focal_core_iqr"
        if len(levels) < 3:
            # quartiles collapsed (few boxes): keep the median, widen to x0.5 and x1.5 of it
            levels = sorted({round(sign * max(0.5 * q2, floor), 2), round(sign * max(q2, floor), 2), round(sign * max(1.5 * q2, floor), 2)})
            source = "fastmri_plus_focal_core_median_widened"
        contrast_by_seq[seq] = levels
        contrast_source[seq] = source
        print(f"    -> bank contrast levels {levels} (sign {sign:+.0f}, floor {floor}, {source})")
    if args.contrast_levels:
        for item in args.contrast_levels.split(";"):
            seq, levels = item.split("=")
            contrast_by_seq[seq.strip()] = sorted(float(x) for x in levels.split(","))
            contrast_source[seq.strip()] = "forced_amendment"
        print(f"  forced contrast levels: { {k: v for k, v in contrast_by_seq.items() if contrast_source.get(k) == 'forced_amendment'} }")
    if args.default_contrasts:
        defaults = sorted(float(x) for x in args.default_contrasts.split(","))
        for stem in file_index:
            seq = sequence_of(stem)
            if seq not in contrast_by_seq:
                sign = signs.get(seq, signs.get(seq.replace("AX", ""), [1]))
                sign = float(sign[0]) if isinstance(sign, (list, tuple)) else float(sign)
                contrast_by_seq[seq] = sorted(round(sign * d, 2) for d in defaults)
                contrast_source[seq] = "default_insufficient_focal_boxes" if len(stats) and (stats["sequence"] == seq).any() else "default_no_annotations"
        print(f"  default contrasts {defaults} applied to sequences without statistics: "
              f"{sorted(k for k, v in contrast_source.items() if v == 'default')}")

    excluded = annotated_slices(annotations)
    manifest = build_manifest(file_index, excluded, bank_cfg, contrast_by_seq, seed=args.seed,
                              middle_slices=args.middle_slices, maps_cache=args.maps_cache)
    if len(manifest):
        manifest["contrast_source"] = manifest["sequence"].map(contrast_source)
    os.makedirs(os.path.dirname(args.manifest_out) or ".", exist_ok=True)
    manifest.to_csv(args.manifest_out, index=False)
    print(f"bank manifest: {len(manifest)} planned insertions -> {args.manifest_out}")


if __name__ == "__main__":
    main()
