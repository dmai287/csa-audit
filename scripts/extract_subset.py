#!/usr/bin/env python3
"""Stream a fastMRI .tar.xz once and extract only wanted volumes, stopping early.

The archives are single-stream xz, so tar cannot seek; listing or extracting
anything costs a full decompression pass. This script pipes `xz -dc` through
Python's streaming tarfile reader, extracts members whose stem is in the wanted
list (optionally balanced by sequence), and terminates the pipe as soon as the
quota is met. It also logs every wanted stem it passes, so later runs know
where files sit in the stream.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tarfile
import time
from collections import Counter

SEQ_RE = re.compile(r"file_brain_(AX[A-Z0-9]+)_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--wanted", required=True, help="CSV/txt with one file stem per line (fastMRI+ brain_file_list.csv)")
    ap.add_argument("--out", required=True, help="output directory (members keep their archive-relative path)")
    ap.add_argument("--quota", default="AXFLAIR=15,AXT1=10,AXT1POST=5",
                    help="per-sequence maximum, e.g. AXFLAIR=15,AXT1=10; sequences not listed are skipped")
    ap.add_argument("--seen-log", default=None, help="append every wanted stem encountered (name, position) here")
    args = ap.parse_args()

    wanted = set()
    with open(args.wanted) as f:
        for line in f:
            stem = line.strip().split(",")[0]
            if stem.startswith("file_brain_"):
                wanted.add(stem)
    quota = {k: int(v) for k, v in (kv.split("=") for kv in args.quota.split(","))}
    got = Counter()
    os.makedirs(args.out, exist_ok=True)
    seen = open(args.seen_log, "a") if args.seen_log else None

    proc = subprocess.Popen(["xz", "-dc", "-T0", args.archive], stdout=subprocess.PIPE, bufsize=1 << 22)
    t0 = time.time(); n_members = 0; n_bytes = 0
    try:
        with tarfile.open(fileobj=proc.stdout, mode="r|") as tar:
            for m in tar:
                n_members += 1; n_bytes += m.size
                if not m.isfile():
                    continue
                stem = os.path.splitext(os.path.basename(m.name))[0]
                if stem not in wanted:
                    continue
                seq = (SEQ_RE.match(stem) or [None, "?"])[1]
                if seen:
                    seen.write(f"{os.path.basename(args.archive)},{stem},{seq},{n_members}\n"); seen.flush()
                if seq not in quota or got[seq] >= quota[seq]:
                    continue
                dst = os.path.join(args.out, m.name)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                src = tar.extractfile(m)
                with open(dst, "wb") as out:
                    for chunk in iter(lambda: src.read(1 << 22), b""):
                        out.write(chunk)
                if os.path.getsize(dst) != m.size:
                    raise RuntimeError(f"size mismatch extracting {m.name}")
                got[seq] += 1
                print(f"[{time.time()-t0:6.0f}s] {stem} ({m.size/1e6:.0f} MB) -> {dict(got)}; streamed {n_bytes/1e9:.1f} GB", flush=True)
                if all(got[s] >= q for s, q in quota.items()):
                    print("quota met; stopping stream", flush=True)
                    break
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        if seen:
            seen.close()
    print(f"done: {sum(got.values())} files, {n_members} members scanned, {n_bytes/1e9:.1f} GB streamed, {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
