#!/usr/bin/env python3
"""Experiment 1: counterfactual lesion stress test over the bank manifest.

Skeleton: iterates the manifest, builds the operator per volume (ESPIRiT maps
from the calibration region), inserts each lesion by Algorithm 1, runs the
configured reconstructors and writes one row per pair with provenance columns.
Requires fastMRI data, the manifest from 02_bank.py and configured weights.
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    except Exception:
        return "unknown"


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="outputs/bank_manifest.csv")
    ap.add_argument("--data", required=True)
    ap.add_argument("--models", default="configs/models.yaml")
    ap.add_argument("--thresholds", default="configs/thresholds.yaml")
    ap.add_argument("--out", default="outputs/exp1.parquet")
    args = ap.parse_args()
    raise SystemExit(
        "Not runnable yet: needs data, manifest and weights. Provenance helpers "
        f"(git {git_commit()[:8]}) are in place; see scripts/00_synthetic_smoke.py for the loop shape.")


if __name__ == "__main__":
    main()
