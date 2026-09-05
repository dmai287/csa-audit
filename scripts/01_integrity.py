#!/usr/bin/env python3
"""Integrity check of fastMRI brain multi-coil volumes before any audit.

Flags zero-filled k-space regions (partial Fourier or asymmetric echo), records
coil counts, matrix sizes and pixel spacing, and writes a CSV manifest. Volumes
with any zero-filled fraction are excluded from the audit.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.io.fastmri import pixel_spacing_mm, read_volume, zero_padding_report  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="directory of fastMRI brain multicoil .h5 files")
    ap.add_argument("--out", default="outputs/integrity.csv")
    args = ap.parse_args()
    rows = []
    for path in sorted(glob.glob(os.path.join(args.data, "*.h5"))):
        vol = read_volume(path)
        rep = zero_padding_report(vol["kspace"])
        try:
            spacing = pixel_spacing_mm(vol["header"])
        except Exception:  # header variants
            spacing = {}
        rows.append({"file": os.path.basename(path), **rep,
                     "spacing_x_mm": spacing.get("x"), "spacing_y_mm": spacing.get("y"),
                     "spacing_z_mm": spacing.get("z"), "acquisition": vol["attrs"].get("acquisition"),
                     "excluded": rep["zero_line_fraction"] > 0 or rep["zero_row_fraction"] > 0})
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    print(df.describe(include="all").T if len(df) else "no files found")
    print(f"{int(df['excluded'].sum()) if len(df) else 0} of {len(df)} volumes excluded")


if __name__ == "__main__":
    main()
