"""Cluster bootstrap, equivalence testing and sample-size rules (paper Section 4.8)."""
from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm

from csa.observer.cho import auc_mann_whitney


def cluster_bootstrap(df: pd.DataFrame, stat: Callable[[pd.DataFrame], float], cluster: str = "subject",
                      n_boot: int = 2000, seed: int = 0, alpha: float = 0.05) -> Tuple[float, float, float]:
    """Resample clusters with replacement; return (estimate, lo, hi) percentile interval."""
    rng = np.random.default_rng(seed)
    ids = df[cluster].unique()
    groups = {k: g for k, g in df.groupby(cluster)}
    est = float(stat(df))
    boots = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(ids, size=len(ids), replace=True)
        boots[b] = stat(pd.concat([groups[k] for k in pick], ignore_index=True))
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    return est, float(lo), float(hi)


def rate_ci(df: pd.DataFrame, indicator: str, **kw) -> Tuple[float, float, float]:
    """Cluster-bootstrap interval for a proportion such as the silent erasure rate."""
    return cluster_bootstrap(df, lambda d: d[indicator].mean(), **kw)


def standardised_mean_difference(x: np.ndarray, y: np.ndarray) -> float:
    sp = np.sqrt(0.5 * (np.var(x, ddof=1) + np.var(y, ddof=1)))
    return float((np.mean(x) - np.mean(y)) / sp)


def tost_equivalence(df: pd.DataFrame, value: str, group: str, margin: float = 0.2,
                     alpha: float = 0.05, **kw) -> dict:
    """Equivalence of a standardised mean difference by the two-one-sided-tests logic:
    declared when the (1 - 2 alpha) cluster-bootstrap interval lies within +/- margin."""
    def smd(d):
        g = d[group].astype(bool)
        return standardised_mean_difference(d.loc[g, value].values, d.loc[~g, value].values)
    est, lo, hi = cluster_bootstrap(df, smd, alpha=2 * alpha, **kw)
    return {"smd": est, "ci_lo": lo, "ci_hi": hi, "margin": margin, "equivalent": bool(lo > -margin and hi < margin)}


def auc_ci(df: pd.DataFrame, score: str, label: str, **kw) -> Tuple[float, float, float]:
    """Cluster-bootstrap interval for the AUC of `score` as a classifier of `label`."""
    def auc(d):
        g = d[label].astype(bool)
        return auc_mann_whitney(d.loc[g, score].values, d.loc[~g, score].values)
    return cluster_bootstrap(df, auc, **kw)


def design_effect(m: float, rho: float) -> float:
    """1 + (m - 1) rho for m pairs per subject and intra-subject correlation rho."""
    return 1.0 + (m - 1.0) * rho


def pairs_per_cell(pi0: float, half_width: float, m: float, rho: float, level: float = 0.95) -> int:
    """Pairs needed so the Wald half-width for a rate pi0 is <= half_width after the design effect."""
    z = norm.ppf(1 - (1 - level) / 2)
    n0 = z ** 2 * pi0 * (1 - pi0) / half_width ** 2
    return int(np.ceil(n0 * design_effect(m, rho)))
