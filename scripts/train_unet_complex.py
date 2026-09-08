#!/usr/bin/env python3
"""Train the two-channel complex U-Net baseline on fastMRI brain data.

Not runnable yet on this machine: needs torch, fastmri and a CUDA device
(see docs/MODELS.md). No lesion-bearing slice may enter training or
validation; the exclusion list comes from fastMRI+ annotations.
"""
from __future__ import annotations

import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-dir", required=True, help="extracted multicoil_train directory")
    ap.add_argument("--annotations", required=True, help="fastMRI+ brain.csv, to exclude annotated slices")
    ap.add_argument("--out", default="checkpoints/unet_complex_brain.pt")
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--chans", type=int, default=32)
    ap.add_argument("--num-pool-layers", type=int, default=4)
    args = ap.parse_args()
    raise SystemExit(
        "Not implemented: needs torch + fastmri + a CUDA device, none present on this "
        "machine (see docs/MODELS.md). Loop shape: for each fully sampled train volume "
        "not in the fastMRI+ exclusion list, apply a random Cartesian mask (configs/masks.yaml), "
        "zero-fill, normalise, and regress the two-channel complex U-Net (in_chans=out_chans=2) "
        "against the fully sampled sensitivity-combined image with an L1 loss.")


if __name__ == "__main__":
    main()
