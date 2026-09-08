#!/usr/bin/env python3
"""Simulation check of the statistical machinery in csa.stats (CPU, synthetic).

1. Coverage of the subject-cluster bootstrap interval for a rate with
   intra-subject correlation, against the naive binomial interval that ignores
   clustering (expected to under-cover).
2. Behaviour of the TOST equivalence rule at true SMD 0 (should declare
   equivalence often at adequate n) and at true SMD 0.5 (should not).
3. The pairs-per-cell rule against the empirical half-width it produces.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.stats.bootstrap import pairs_per_cell, rate_ci, tost_equivalence  # noqa: E402


def clustered_binary(n_subjects, m, pi, rho, rng):
    """Beta-binomial subjects: E[p]=pi, ICC=rho."""
    kappa = (1 - rho) / rho
    p = rng.beta(pi * kappa, (1 - pi) * kappa, size=n_subjects)
    subj = np.repeat(np.arange(n_subjects), m)
    y = rng.random(n_subjects * m) < np.repeat(p, m)
    return pd.DataFrame({"subject": subj, "silent_erased": y.astype(int)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--n-boot", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    pi, rho, m = 0.05, 0.10, 12
    n_subj = int(np.ceil(pairs_per_cell(pi, 0.02, m, rho) / m))
    print(f"design: pi0={pi}, rho={rho}, m={m} -> pairs_per_cell={pairs_per_cell(pi, 0.02, m, rho)}, subjects={n_subj}")
    cover_cb, cover_naive, widths = 0, 0, []
    for k in range(args.reps):
        df = clustered_binary(n_subj, m, pi, rho, rng)
        est, lo, hi = rate_ci(df, "silent_erased", n_boot=args.n_boot, seed=k)
        cover_cb += lo <= pi <= hi
        widths.append((hi - lo) / 2)
        n = len(df); se = np.sqrt(est * (1 - est) / n); z = norm.ppf(0.975)
        cover_naive += est - z * se <= pi <= est + z * se
    print(f"cluster-bootstrap 95% CI coverage: {cover_cb/args.reps:.2f}   (naive binomial: {cover_naive/args.reps:.2f})")
    print(f"empirical half-width: mean {np.mean(widths):.4f} (target <= 0.02)")

    for smd_true in (0.0, 0.5):
        declared = 0
        for k in range(args.reps // 2):
            n_s = 40
            subj = np.repeat(np.arange(n_s), 10)
            grp = rng.random(len(subj)) < 0.5
            u = np.repeat(rng.normal(0, 0.3, n_s), 10)
            val = u + rng.normal(0, 1, len(subj)) + smd_true * grp
            df = pd.DataFrame({"subject": subj, "erased": grp.astype(int), "psnr": val})
            declared += tost_equivalence(df, "psnr", "erased", margin=0.2, n_boot=args.n_boot, seed=k)["equivalent"]
        print(f"TOST (margin 0.2): true SMD {smd_true} -> equivalence declared in {declared/(args.reps//2):.2f} of runs")


if __name__ == "__main__":
    main()
