"""One counterfactual pair -> one audit row (Algorithm 1 with all readouts)."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable, Dict, Optional

import numpy as np

from csa.lesion.insert import counterfactual_kspace
from csa.metrics.image import annulus, psnr, roi_cnr, roi_ssim, ssim
from csa.nullspace.projectors import measurement_fraction, measurement_gain, transfer_statistics
from csa.observer.cho import CHO, extract_roi, per_lesion_z
from csa.physics.operator import SenseOperator


@dataclass
class PairSpec:
    subject: str
    slice_index: int
    site_row: int
    site_col: int
    volume_mm3: float
    contrast: float
    model: str
    acceleration: int
    mask_type: str
    seed: int = 0


Reconstructor = Callable[[np.ndarray, np.ndarray, np.ndarray, int], np.ndarray]


def audit_pair(spec: PairSpec, op: SenseOperator, y_full: np.ndarray, X: np.ndarray,
               ell: np.ndarray, reconstruct: Reconstructor, lam: float,
               roi_size: int, cho: Optional[CHO] = None, t0_recon: Optional[np.ndarray] = None,
               t0_ref: Optional[np.ndarray] = None, z_ref: Optional[float] = None,
               radius_px: float = 2.0, cg_tol: float = 1e-8,
               cho_ref: Optional[CHO] = None) -> Dict:
    """Run both arms and return one row of readouts.

    `reconstruct(y_masked, mask, maps, seed)` must return a complex image.
    `cho` is an observer fitted on reconstructed patches and `t0_recon` its
    lesion-absent statistics; `cho_ref` (default `cho`) is the observer for the
    fully sampled reference with lesion-absent statistics `t0_ref`. The same
    channels are used for both; the templates are fitted per domain.
    """
    y = op.undersample(y_full)
    y_cf = counterfactual_kspace(op, y, ell)
    Xh = reconstruct(y, op.mask, op.maps, spec.seed)
    Xh_cf = reconstruct(y_cf, op.mask, op.maps, spec.seed)
    dX = Xh_cf - Xh
    X_cf = X + ell
    center = (spec.site_row, spec.site_col)
    disc = np.hypot(*np.mgrid[0:X.shape[0], 0:X.shape[1]] - np.array(center)[:, None, None]) <= radius_px
    ring = annulus(X.shape, center, radius_px + 2, radius_px + 8)

    t_R, t_N = transfer_statistics(op, dX, ell, lam, tol=cg_tol)
    row = asdict(spec)
    row.update({
        "residual_factual": op.relative_residual(Xh, y),
        "residual_cf": op.relative_residual(Xh_cf, y_cf),
        "measurement_gain": measurement_gain(op, ell),
        "mu_lambda": measurement_fraction(op, ell, lam, tol=cg_tol),
        "t_R": t_R,
        "t_N": t_N,
        "psnr_factual": psnr(X, Xh),
        "psnr_cf": psnr(X_cf, Xh_cf),
        "ssim_factual": ssim(X, Xh),
        "ssim_cf": ssim(X_cf, Xh_cf),
        "roi_cnr_cf": roi_cnr(Xh_cf, disc, ring),
        "roi_cnr_ref": roi_cnr(X_cf, disc, ring),
        "roi_ssim_cf": roi_ssim(X_cf, Xh_cf, center, roi_size),
        "lesion_energy": float(np.vdot(ell, ell).real),
    })
    if cho is not None and t0_recon is not None:
        t_recon = float(cho.decide(extract_roi(Xh_cf, center, roi_size)[None])[0])
        row["z_recon"] = per_lesion_z(t_recon, t0_recon)
        if t0_ref is not None:
            t_ref = float((cho_ref or cho).decide(extract_roi(X_cf, center, roi_size)[None])[0])
            row["z_ref"] = per_lesion_z(t_ref, t0_ref)
        elif z_ref is not None:
            row["z_ref"] = float(z_ref)
    return row
