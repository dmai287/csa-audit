#!/usr/bin/env python3
"""Figures from the per-lesion tables. Reads only CSV outputs, never intermediate arrays.

validation mode (default): figures from the model-independent characterization
and the classical-reconstructor demonstrations. results mode: reserved for the
confirmatory Experiment 1 table (learned models), not implemented until it exists.
"""
from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VOL_ORDER = [1.0, 5.0, 10.0, 27.0]


def fig_measurement(ch, out):
    """kappa^2 and mu_lambda against acceleration, one line per lesion volume, with 1/R."""
    u = ch.drop_duplicates(["file", "slice", "site_row", "site_col", "volume_mm3", "acceleration"])
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharex=True)
    Rs = sorted(u.acceleration.unique())
    for ax, col, name in zip(axes, ["kappa2", "mu_lambda"], [r"measurement gain $\kappa^2(\ell)$", r"measured fraction $\mu_\lambda(\ell)$"]):
        for v in VOL_ORDER:
            g = u[u.volume_mm3 == v].groupby("acceleration")[col]
            m, lo, hi = g.median(), g.quantile(0.25), g.quantile(0.75)
            ax.errorbar(m.index, m.values, yerr=[m.values - lo.values, hi.values - m.values], marker="o", capsize=3, label=f"{v:g} mm$^3$")
        ax.plot(Rs, [1 / r for r in Rs], "k--", lw=1, label=r"single-coil $1/R$")
        ax.set_xlabel("acceleration $R$"); ax.set_ylabel(name); ax.set_xticks(Rs); ax.set_ylim(0, 1); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, ncol=2)
    fig.suptitle("Model-independent: how much of each planned lesion the acquisition measures (median, IQR over lesions)", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def fig_zref(ch, z_det, out):
    """Reference detectability by volume and contrast, with the pre-registered z_det."""
    zr = ch[ch.acceleration == ch.acceleration.min()].dropna(subset=["z_ref"])
    seqs = sorted(zr.sequence.unique())
    fig, axes = plt.subplots(1, len(seqs), figsize=(4.5 * len(seqs), 3.6), squeeze=False)
    for ax, seq in zip(axes[0], seqs):
        g = zr[zr.sequence == seq]
        for c in sorted(g.contrast.unique()):
            gg = g[g.contrast == c].groupby("volume_mm3")["z_ref"]
            ax.errorbar(gg.median().index, gg.median().values, yerr=[gg.median().values - gg.quantile(0.25).values, gg.quantile(0.75).values - gg.median().values],
                        marker="o", capsize=3, label=f"contrast {c:+.2f}")
        ax.axhline(z_det, color="k", ls="--", lw=1, label=f"$z_{{det}}$ = {z_det}")
        ax.set_xscale("log"); ax.set_xlabel("lesion volume (mm$^3$)"); ax.set_ylabel("reference detectability $z_{ref}$"); ax.set_title(seq, fontsize=10)
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle("Detectability of the inserted lesion in the fully sampled reference (cross-fitted CHO; median, IQR)", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def fig_exp2(e2, out):
    """Residual insensitivity: null-space edits vs the residual, and error split."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    ax = axes[0]
    ax.scatter(e2.edit_A_norm_rel, e2.residual_change_from_edit, s=18)
    lim = [min(e2.residual_change_from_edit.min(), e2.edit_A_norm_rel.min()) * 0.5, e2.edit_A_norm_rel.max() * 2]
    ax.plot(lim, lim, "k--", lw=1, label="Proposition 1 bound")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"$\|A\delta\|_2 / \|Y\|_2$ (null-space lesion edit)"); ax.set_ylabel("change in relative residual")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[1]
    w = 0.25; x = np.arange(len(e2))
    ax.bar(x - w, e2.residual_full, w, label="reconstruction"); ax.bar(x, e2.residual_without_null_error, w, label="minus null-space error")
    ax.bar(x + w, e2.residual_without_range_error, w, label="minus measured error")
    ax.set_xlabel("lesion case"); ax.set_ylabel("relative k-space residual"); ax.legend(fontsize=7); ax.grid(alpha=0.3, axis="y")
    fig.suptitle("Residual insensitivity on real acquisitions (classical reconstructor; validation)", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def fig_exp3(e3, out):
    """mu_lambda under acquisition shift, per condition, relative to the equispaced baseline."""
    base = e3[e3.condition == "equispaced_offset_0"].set_index(["file", "slice", "site_row", "site_col", "volume_mm3"])["mu_lambda"]
    d = e3.set_index(["file", "slice", "site_row", "site_col", "volume_mm3"])
    d["rel"] = d["mu_lambda"] / base.reindex(d.index).values
    order = [c for c in ["equispaced_offset_1", "equispaced_offset_2", "equispaced_offset_3", "random_seed_0", "random_seed_1", "random_seed_2",
                         "noise_x2", "noise_x4", "coils_first_half", "coils_second_half", "coils_alternate"] if c in set(d.condition)]
    fig, ax = plt.subplots(figsize=(9, 3.6))
    data = [d[d.condition == c]["rel"].dropna().values for c in order]
    ax.boxplot(data, labels=[c.replace("_", "\n") for c in order]); ax.axhline(1, color="k", ls="--", lw=1)
    ax.set_ylabel(r"$\mu_\lambda$ relative to equispaced baseline"); ax.tick_params(axis="x", labelsize=7); ax.grid(alpha=0.3, axis="y")
    fig.suptitle("Measurement-side acquisition shift (model-independent; validation)", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def fig_explore(ex, z_det, out):
    """z_ref over the extended volume x contrast grid, thick slice vs thin-slice equivalent."""
    seqs = sorted(ex.sequence.unique())
    fig, axes = plt.subplots(2, len(seqs), figsize=(4.6 * len(seqs), 6.4), squeeze=False)
    for j, seq in enumerate(seqs):
        for i, (fill, name) in enumerate(((True, "thick slice (as acquired)"), (False, "thin-slice equivalent"))):
            ax = axes[i][j]; g = ex[(ex.sequence == seq) & (ex.with_fill == fill)]
            for c in sorted(g.contrast.unique(), key=abs):
                gg = g[g.contrast == c].groupby("volume_mm3")["z_median"].mean()
                ax.plot(gg.index, gg.values, marker="o", label=f"contrast {c:+.2f}")
            ax.axhline(z_det, color="k", ls="--", lw=1)
            ax.set_xscale("log"); ax.set_title(f"{seq}: {name}", fontsize=9); ax.grid(alpha=0.3)
            ax.set_xlabel("lesion volume (mm$^3$)"); ax.set_ylabel("median $z_{ref}$"); ax.legend(fontsize=7)
    fig.suptitle(f"Reference detectability over an extended grid (dashed: pre-registered $z_{{det}}$ = {z_det})", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def fig_dryrun(df, z_det, z_miss, out):
    """Classical dry run: detectability reference vs reconstruction, and PSNR change per pair."""
    models = sorted(df.model.unique())
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    ax = axes[0]
    for m in models:
        g = df[df.model == m].dropna(subset=["z_ref", "z_recon"])
        ax.scatter(g.z_ref, g.z_recon, s=10, alpha=0.6, label=m)
    lim = [min(df.z_ref.min(), df.z_recon.min()), max(df.z_ref.max(), df.z_recon.max())]
    ax.plot(lim, lim, "k--", lw=1); ax.axvline(z_det, color="k", ls=":", lw=1); ax.axhline(z_miss, color="k", ls=":", lw=1)
    ax.set_xlabel("$z_{ref}$ (fully sampled reference)"); ax.set_ylabel("$z_{recon}$ (reconstruction)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("per-lesion detectability (dotted: $z_{det}$, $z_{miss}$)", fontsize=9)
    ax = axes[1]
    for m in models:
        g = df[df.model == m]
        ax.hist(g.psnr_cf - g.psnr_factual, bins=40, alpha=0.6, label=m)
    ax.set_xlabel("PSNR(lesion present) - PSNR(lesion absent), dB"); ax.set_ylabel("pairs"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("what the global score sees of the lesion", fontsize=9)
    ax = axes[2]
    for m in models:
        g = df[df.model == m]
        ax.scatter(g.mu_lambda, g.t_N, s=10, alpha=0.6, label=m)
    ax.set_xlabel(r"measured fraction $\mu_\lambda$"); ax.set_ylabel("null-space transfer $t_N$"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("what the reconstruction supplied beyond the measurement", fontsize=9)
    fig.suptitle("Experiment 1 machinery, classical reconstructors at R=8 (dry run; not a confirmatory result)", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default="outputs/annotated_val0")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--thresholds", default="configs/thresholds.yaml")
    args = ap.parse_args()
    out = args.out_dir or os.path.join(args.run_dir, "figures"); os.makedirs(out, exist_ok=True)
    import yaml
    z_det = yaml.safe_load(open(args.thresholds)).get("z_det") or 3.0
    ch = pd.read_csv(os.path.join(args.run_dir, "bank_characterization.csv"))
    fig_measurement(ch, os.path.join(out, "val_measurement_gain_fraction.png"))
    fig_zref(ch, z_det, os.path.join(out, "val_reference_detectability.png"))
    p2 = os.path.join(args.run_dir, "exp2_physics_demo.csv")
    if os.path.exists(p2):
        fig_exp2(pd.read_csv(p2), os.path.join(out, "val_exp2_residual_insensitivity.png"))
    p3 = os.path.join(args.run_dir, "exp3_measurement_side.csv")
    if os.path.exists(p3):
        fig_exp3(pd.read_csv(p3), os.path.join(out, "val_exp3_acquisition_shift.png"))
    pe = os.path.join(args.run_dir, "explore_reference_detectability.csv")
    if os.path.exists(pe):
        fig_explore(pd.read_csv(pe), z_det, os.path.join(out, "val_explore_reference_detectability.png"))
    pd1 = os.path.join(args.run_dir, "exp1_dryrun_classical.csv")
    if os.path.exists(pd1):
        thr = yaml.safe_load(open(args.thresholds))
        fig_dryrun(pd.read_csv(pd1), thr.get("z_det") or 3.0, thr.get("z_miss") or 2.0, os.path.join(out, "dryrun_exp1_classical.png"))
    print("figures written to", out, ":", sorted(os.listdir(out)))


if __name__ == "__main__":
    main()
