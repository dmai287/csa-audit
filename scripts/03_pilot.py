#!/usr/bin/env python3
"""Pilot analysis: intra-subject correlation and pairs per cell (paper Section 4.8).

Reads a pilot audit table (CSV with columns subject, silent_erased) and returns
the number of pairs per model-by-R cell needed for the pre-registered interval
half-width, using the design effect 1 + (m - 1) rho.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.stats.bootstrap import design_effect, pairs_per_cell  # noqa: E402


def icc_binary(df: pd.DataFrame, value: str, cluster: str) -> float:
    """ANOVA-type intraclass correlation for a binary indicator."""
    g = df.groupby(cluster)[value]
    k = g.size().mean()
    msb = g.mean().var(ddof=1) * k
    msw = (g.var(ddof=1).fillna(0) * (g.size() - 1)).sum() / max((g.size() - 1).sum(), 1)
    return float(max((msb - msw) / (msb + (k - 1) * msw), 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", required=True, help="pilot audit table (CSV)")
    ap.add_argument("--pi0", type=float, required=True)
    ap.add_argument("--half-width", type=float, default=0.02)
    args = ap.parse_args()
    df = pd.read_csv(args.pilot)
    rho = icc_binary(df, "silent_erased", "subject")
    m = df.groupby("subject").size().mean()
    n = pairs_per_cell(args.pi0, args.half_width, m, rho)
    print(f"rho = {rho:.3f}, m = {m:.1f}, design effect = {design_effect(m, rho):.2f}, pairs per cell = {n}")


if __name__ == "__main__":
    main()
