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


def appearance_stats(annotations: pd.DataFrame, file_index: dict, ring_px: int = 8) -> pd.DataFrame:
    """For each annotated bounding box present locally: RSS magnitude in the box
    vs. an annulus around it, and the resulting contrast ratio."""
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
            rows.append({"file": r["file"], "slice": sl, "sequence": r["sequence"],
                        "label": r["label"], "width": w, "height": h,
                        "box_mag": box_mag, "background_mag": bg_mag,
                        "contrast_ratio": (box_mag - bg_mag) / bg_mag if bg_mag > 0 else np.nan})
    finally:
        for f in open_files.values():
            f.close()
    cols = ["file", "slice", "sequence", "label", "width", "height", "box_mag", "background_mag", "contrast_ratio"]
    return pd.DataFrame(rows, columns=cols)


def annotated_slices(annotations: pd.DataFrame) -> set:
    """(file, slice) pairs to exclude from insertion, regardless of sequence."""
    return set(zip(annotations["file"], annotations["slice"].astype(int)))


def build_manifest(file_index: dict, excluded: set, bank_cfg: dict, contrast_by_seq: dict,
                   seed: int = 0, middle_slices: int = 0) -> pd.DataFrame:
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
            candidates = np.argwhere((a > 0.25) & (a < 0.9))
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
    for seq, g in (stats.groupby("sequence") if len(stats) else []):
        q1, q3 = g["contrast_ratio"].quantile([0.25, 0.75])
        contrast_by_seq[seq] = sorted({round(q1, 2), round((q1 + q3) / 2, 2), round(q3, 2)})
        contrast_source[seq] = "fastmri_plus_iqr"
        print(f"  {seq}: n={len(g)}, contrast IQR [{q1:.2f}, {q3:.2f}]")
    if args.default_contrasts:
        defaults = sorted(float(x) for x in args.default_contrasts.split(","))
        for stem in file_index:
            seq = sequence_of(stem)
            if seq not in contrast_by_seq:
                contrast_by_seq[seq] = defaults
                contrast_source[seq] = "default"
        print(f"  default contrasts {defaults} applied to sequences without statistics: "
              f"{sorted(k for k, v in contrast_source.items() if v == 'default')}")

    excluded = annotated_slices(annotations)
    manifest = build_manifest(file_index, excluded, bank_cfg, contrast_by_seq, seed=args.seed,
                              middle_slices=args.middle_slices)
    if len(manifest):
        manifest["contrast_source"] = manifest["sequence"].map(contrast_source)
    os.makedirs(os.path.dirname(args.manifest_out) or ".", exist_ok=True)
    manifest.to_csv(args.manifest_out, index=False)
    print(f"bank manifest: {len(manifest)} planned insertions -> {args.manifest_out}")


if __name__ == "__main__":
    main()
