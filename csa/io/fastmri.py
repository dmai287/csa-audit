"""fastMRI brain multi-coil file access and integrity checks.

Files are HDF5 with `kspace` (slices, coils, H, W) complex64 and an
`ismrmrd_header` XML string carrying field of view and matrix size.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import Dict

import h5py
import numpy as np


def read_volume(path: str) -> Dict:
    with h5py.File(path, "r") as f:
        ks = np.asarray(f["kspace"])
        header = f["ismrmrd_header"][()]
        attrs = dict(f.attrs)
    header = header.decode() if isinstance(header, bytes) else str(header)
    return {"kspace": ks, "header": header, "attrs": attrs}


def pixel_spacing_mm(header_xml: str) -> Dict[str, float]:
    """Best-effort parse of encodedSpace field of view and matrix size."""
    root = ET.fromstring(header_xml)
    ns = {"i": root.tag.split("}")[0].strip("{")} if root.tag.startswith("{") else {}
    def find(path):
        node = root.find(path, ns) if ns else root.find(path.replace("i:", ""))
        return None if node is None else float(node.text)
    fov = {k: find(f".//i:encoding/i:encodedSpace/i:fieldOfView_mm/i:{k}") for k in ("x", "y", "z")}
    mat = {k: find(f".//i:encoding/i:encodedSpace/i:matrixSize/i:{k}") for k in ("x", "y", "z")}
    out = {}
    for k in ("x", "y", "z"):
        if fov[k] and mat[k]:
            out[k] = fov[k] / mat[k]
    return out


def zero_padding_report(kspace: np.ndarray, tol: float = 0.0) -> Dict:
    """Fraction of phase-encode lines and readout rows that are exactly zero in every coil.

    Non-zero fractions indicate partial Fourier or zero-filled regions, in which case
    the 'fully sampled' reference was not fully measured and the volume is excluded.
    """
    a = np.abs(kspace)
    zero_lines = (a.max(axis=(0, 1, 2)) <= tol).mean()
    zero_rows = (a.max(axis=(0, 1, 3)) <= tol).mean()
    return {"zero_line_fraction": float(zero_lines), "zero_row_fraction": float(zero_rows),
            "n_slices": int(kspace.shape[0]), "n_coils": int(kspace.shape[1]),
            "shape": tuple(int(s) for s in kspace.shape[2:])}
