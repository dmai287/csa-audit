#!/usr/bin/env python3
"""Train the two-channel (real, imaginary) image-domain U-Net baseline.

Input: the SENSE-combined zero-filled image of a randomly masked slice (two
channels, normalised by its own 99th-percentile magnitude). Target: the
fully sampled sensitivity-combined reference in the same normalisation. Loss:
L1 on the two channels. No lesion-bearing slice enters training: slices with
any fastMRI+ annotation are excluded, as are all validation hosts.

Coil maps come from csa.physics.espirit (sigpy) and are cached under
--maps-cache with the same naming as the audit, so the model sees exactly the
combination the audit will use at inference.

GPU: torch with CUDA. CPU: works for a smoke test (--limit-slices, --epochs 1).
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import time

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from csa.physics.espirit import espirit_maps  # noqa: E402
from csa.physics.masks import default_center_fraction, equispaced_mask, random_mask  # noqa: E402
from csa.physics.operator import SenseOperator, ifft2c, sense_combine  # noqa: E402


def load_maps(path, sl, cache_dir):
    key = f"{os.path.splitext(os.path.basename(path))[0]}_s{sl:02d}.npz"
    cp = os.path.join(cache_dir, key)
    if os.path.exists(cp):
        return np.load(cp)["maps"].astype(np.complex128)
    with h5py.File(path, "r") as f:
        ks = np.asarray(f["kspace"][sl])
    maps = espirit_maps(ks, calib_width=24)
    os.makedirs(cache_dir, exist_ok=True)
    np.savez_compressed(cp, maps=maps.astype(np.complex64))
    return maps.astype(np.complex128)


def make_example(path, sl, cache_dir, rng, accelerations=(4, 6, 8)):
    with h5py.File(path, "r") as f:
        ks = np.asarray(f["kspace"][sl]).astype(np.complex128)
    maps = load_maps(path, sl, cache_dir)
    X = sense_combine(ifft2c(ks), maps)
    R = int(rng.choice(accelerations))
    cf = default_center_fraction(R)
    mask = equispaced_mask(ks.shape[-1], R, cf, offset=int(rng.integers(R))) if rng.random() < 0.5 else random_mask(ks.shape[-1], R, cf, seed=int(rng.integers(1 << 30)))
    op = SenseOperator(maps, mask)
    zf = op.adjoint(op.undersample(ks))
    scale = np.quantile(np.abs(zf), 0.99) or 1.0
    x = np.stack([zf.real, zf.imag], 0) / scale
    y = np.stack([X.real, X.imag], 0) / scale
    return x.astype(np.float32), y.astype(np.float32)


def slice_list(train_dir, annotations, exclude_files, limit=0, seed=0):
    ann = pd.read_csv(annotations) if annotations and os.path.exists(annotations) else None
    excluded = set(zip(ann["file"], ann["slice"].astype(int))) if ann is not None else set()
    items = []
    for dp, _, fs in os.walk(train_dir):
        for f in sorted(fs):
            if not f.endswith(".h5"):
                continue
            stem = os.path.splitext(f)[0]
            if stem in exclude_files:
                continue
            with h5py.File(os.path.join(dp, f), "r") as h:
                n = h["kspace"].shape[0]
            for sl in range(2, n - 2):  # drop the outermost slices
                if (stem, sl) not in excluded:
                    items.append((os.path.join(dp, f), sl))
    random.Random(seed).shuffle(items)
    return items[:limit] if limit else items


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-dir", required=True)
    ap.add_argument("--annotations", default="", help="fastMRI+ brain.csv; annotated slices are excluded")
    ap.add_argument("--exclude-hosts", default="", help="text file of file stems never to train on (validation hosts)")
    ap.add_argument("--maps-cache", default="cache/maps")
    ap.add_argument("--out", default="checkpoints/unet_complex_brain.pt")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--chans", type=int, default=32)
    ap.add_argument("--num-pool-layers", type=int, default=4)
    ap.add_argument("--limit-slices", type=int, default=0)
    ap.add_argument("--val-fraction", type=float, default=0.05)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import torch
    from fastmri.models import Unet
    torch.manual_seed(args.seed)
    device = torch.device(args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu")
    exclude = set(open(args.exclude_hosts).read().split()) if args.exclude_hosts and os.path.exists(args.exclude_hosts) else set()
    items = slice_list(args.train_dir, args.annotations, exclude, args.limit_slices, args.seed)
    n_val = max(1, int(len(items) * args.val_fraction)); val_items, train_items = items[:n_val], items[n_val:]
    print(f"{len(train_items)} training slices, {len(val_items)} validation slices, device {device}", flush=True)
    model = Unet(in_chans=2, out_chans=2, chans=args.chans, num_pool_layers=args.num_pool_layers).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.StepLR(opt, step_size=max(args.epochs // 3, 1), gamma=0.3)
    rng = np.random.default_rng(args.seed)
    best = float("inf")
    for epoch in range(args.epochs):
        model.train(); t0 = time.time(); losses = []
        random.Random(args.seed + epoch).shuffle(train_items)
        for path, sl in train_items:
            x, y = make_example(path, sl, args.maps_cache, rng)
            xt = torch.from_numpy(x)[None].to(device); yt = torch.from_numpy(y)[None].to(device)
            # pad to a multiple of 16 for the pooling layers
            H, W = xt.shape[-2:]; ph, pw = (-H) % 16, (-W) % 16
            xt = torch.nn.functional.pad(xt, (0, pw, 0, ph)); yt = torch.nn.functional.pad(yt, (0, pw, 0, ph))
            loss = torch.nn.functional.l1_loss(model(xt), yt)
            opt.zero_grad(); loss.backward(); opt.step(); losses.append(loss.item())
        sched.step()
        model.eval(); vl = []
        with torch.no_grad():
            for path, sl in val_items:
                x, y = make_example(path, sl, args.maps_cache, rng)
                xt = torch.from_numpy(x)[None].to(device); yt = torch.from_numpy(y)[None].to(device)
                H, W = xt.shape[-2:]; ph, pw = (-H) % 16, (-W) % 16
                xt = torch.nn.functional.pad(xt, (0, pw, 0, ph)); yt = torch.nn.functional.pad(yt, (0, pw, 0, ph))
                vl.append(torch.nn.functional.l1_loss(model(xt), yt).item())
        v = float(np.mean(vl))
        print(f"epoch {epoch+1}/{args.epochs}: train L1 {np.mean(losses):.4f}, val L1 {v:.4f}, {time.time()-t0:.0f}s", flush=True)
        if v < best:
            best = v; os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
            torch.save({"state_dict": model.state_dict(), "chans": args.chans, "num_pool_layers": args.num_pool_layers,
                        "in_chans": 2, "out_chans": 2, "epoch": epoch + 1, "val_l1": v}, args.out)
    print(f"best val L1 {best:.4f} -> {args.out}")


if __name__ == "__main__":
    main()
