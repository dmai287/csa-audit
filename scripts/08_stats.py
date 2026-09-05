#!/usr/bin/env python3
"""Rates with cluster-bootstrap intervals from an audit table (works on the smoke output).

For the real analysis the erasure indicators are derived from the pre-registered
thresholds in configs/thresholds.yaml; here they default to illustrative values
and are labelled as such.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.audit.erasure import is_erased  # noqa: E402
from csa.stats.bootstrap import rate_ci, tost_equivalence  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="outputs/smoke.csv")
    ap.add_argument("--thresholds", default="configs/thresholds.yaml")
    args = ap.parse_args()
    df = pd.read_csv(args.table)
    thr = yaml.safe_load(open(args.thresholds))
    z_det, z_miss = thr.get("z_det"), thr.get("z_miss")
    if z_det is None or z_miss is None:
        z_det, z_miss = 3.0, 1.645
        print("WARNING: thresholds not pre-registered; using illustrative z_det=3.0, z_miss=1.645")
    df["erased"] = [is_erased(a, b, z_det, z_miss) for a, b in zip(df["z_ref"], df["z_recon"])]
    for (model, R), g in df.groupby(["model", "acceleration"]):
        est, lo, hi = rate_ci(g, "erased", n_boot=200)
        print(f"{model:12s} R={R}: erasure rate {est:.3f} [{lo:.3f}, {hi:.3f}] (n={len(g)})")
    if df["erased"].nunique() == 2:
        print("H2 (PSNR equivalence between erased and preserved):", tost_equivalence(df, "psnr_cf", "erased", n_boot=200))


if __name__ == "__main__":
    main()
