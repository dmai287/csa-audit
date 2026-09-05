#!/usr/bin/env python3
"""Build the lesion-bank manifest from fastMRI+ annotations and the integrity manifest.

TODO before use: derive contrast levels per sequence from the annotated lesions
(interquartile range of lesion-to-surround intensity ratio), sample insertion
sites by tissue class, and exclude slices with any annotation in the audited
region. Writes outputs/bank_manifest.csv with one row per planned pair.
"""
from __future__ import annotations

import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations", required=True, help="fastMRI+ brain annotation CSV")
    ap.add_argument("--integrity", default="outputs/integrity.csv")
    ap.add_argument("--config", default="configs/bank.yaml")
    ap.add_argument("--out", default="outputs/bank_manifest.csv")
    args = ap.parse_args()
    raise SystemExit("Not implemented yet: needs fastMRI+ CSV in hand. See docstring.")


if __name__ == "__main__":
    main()
