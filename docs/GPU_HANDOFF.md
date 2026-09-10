# Running the confirmatory stage on another machine

Everything CPU-feasible has been run here (see `docs/VALIDATION.md`). What
remains needs a CUDA GPU: training the complex U-Net, the three learned
reconstructors in Experiment 1, the pilot, and Experiments 2 to 4 on the
learned models. This is the hand-off.

## 0. Order of operations (fixed by the pre-registration)

1. Deposit pre-registration v1 and the amendment addendum on OSF
   (`prereg/addendum_v1_DRAFT.md`, filled in). Cut and push a tagged release
   so the addendum cites the exact code (`git tag v0.2.0 && git push --tags`;
   Zenodo archives it if the repository is public at that moment).
2. On the GPU machine: environment, data, checkpoints, U-Net training,
   adapter checks (sections 1 to 5).
3. The 20-subject pilot (`scripts/03_pilot.py`), then post the realised
   sample size as addendum 2.
4. Experiments 1 to 4 with the learned models.

## 1. Environment

```
git clone git@github.com:dmai287/csa-audit.git && cd csa-audit
python -m venv .venv && source .venv/bin/activate      # Python 3.10 or 3.11
pip install -r requirements-gpu.txt && pip install -e .
python -m pytest -q                                      # 20 gate tests on synthetic data
git clone https://github.com/utcsilab/csgm-mri-langevin.git third_party/csgm-mri-langevin
```

## 2. Data (never uploaded from the Mac; pull directly)

The fastMRI links are signed S3 URLs valid until **2026-12-07**. They live
only in `~/.config/fastmri/urls.txt` on the Mac (one `name url` per line) and
must be copied to the GPU machine by hand; they are not in the repository.
`scripts/download_fastmri.sh` in `/Volumes/T9/fastMRI_brain/` is the resumable
downloader; copy it too, or re-use `curl -C -` per line.

Needed for the confirmatory run:
- `brain_multicoil_val_batch_{0,1,2}.tar.xz` (hosts; 296 GB)
- `brain_multicoil_train_batch_0.tar.xz` (U-Net training; 106 GB)
- fastMRI+ annotations: `git clone https://github.com/microsoft/fastmri-plus`

The archives are single-stream xz (no seeking). Extract selectively with
`scripts/extract_subset.py` (streams once, stops at a quota) or fully if disk
allows (val batch 0 alone is 187 GB uncompressed).

Alternatively attach the T9 drive: `/Volumes/T9/fastMRI_brain/extracted/multicoil_val`
holds the 27 host volumes already used (~11 GB) and `cache/maps/` the ESPIRiT
maps for their slices.

## 3. Checkpoints

```
mkdir -p checkpoints
curl -o checkpoints/varnet_brain_leaderboard_state_dict.pt \
  https://dl.fbaipublicfiles.com/fastMRI/trained_models/varnet/brain_leaderboard_state_dict.pt   # 114 MB
pip install gdown && gdown "https://drive.google.com/uc?id=1vAIXf8n67yEAPmH2I9qiDWzmq9fGKPYL" -O csgm.tar.gz
tar -xzf csgm.tar.gz -C checkpoints    # -> checkpoints/ncsnv2-mri-mvue/logs/mri-mvue/checkpoint_100000.pth (1.5 GB)
```
Update the paths in `configs/models.yaml`.

## 4. The complex U-Net (the only model trained here)

```
python scripts/train_unet_complex.py --train-dir data/multicoil_train \
  --annotations fastmri-plus/Annotations/brain.csv --exclude-hosts runs/hosts_val.txt \
  --maps-cache cache/maps --epochs 30 --out checkpoints/unet_complex_brain.pt
```
Roughly one GPU-day on train batch 0. Annotated slices and every validation
host are excluded from training.

## 5. Adapter checks before any experiment

- VarNet: `csa/recon/varnet.py` returns the SENSE-combined complex image from
  the last cascade with the network's own maps. Check on a few slices that
  `abs(adapter output)` tracks `adapter.public_magnitude(...)` (RSS) up to the
  combination rule, and that the residual on sampled lines is small.
- U-Net: `csa/recon/unet.py` reads the architecture from the checkpoint.
- Diffusion: `csa/recon/diffusion.py` wraps the released `LangevinOptimizer`
  unchanged. Construct it with `repo_dir=third_party/csgm-mri-langevin`,
  `checkpoint=...checkpoint_100000.pth`, and `sequence` matching the host
  (selects the released config: T2 -> 384 grid, FLAIR/T1 -> 320 grid).
  First run: reconstruct 3 fully sampled slices at R=4 and check PSNR against
  the reference is in the range the paper reports (Jalal et al 2021, brain);
  time one reconstruction to size the budget. Same seed for both arms of a
  pair (`seed`), replicates on the 10 % subsample only.

## 6. What to run, in order

```
# model-independent characterization for the confirmatory hosts (CPU, resumable)
python scripts/02_bank.py --data data/multicoil_val --middle-slices 6 \
  --contrast-levels "AXFLAIR=0.24,1.0,1.5;AXT1=-0.21,-0.6,-1.0" --default-contrasts 0.3,0.6,1.0 \
  --stats-out outputs/confirm/appearance_stats.csv --manifest-out outputs/confirm/bank_manifest.csv
python scripts/02b_characterize_bank.py --manifest outputs/confirm/bank_manifest.csv --data data/multicoil_val \
  --out outputs/confirm/bank_characterization.csv --workers 16 --cg-maxiter 200
# pilot (20 subjects, rho only), then addendum 2
python scripts/04_exp1.py --manifest outputs/pilot/bank_manifest.csv --characterization outputs/pilot/bank_characterization.csv \
  --data data/pilot_hosts --models e2e_varnet,unet_complex,diffusion_langevin --accelerations 6,8 --label pilot --out-dir outputs/pilot
python scripts/03_pilot.py --pilot outputs/pilot/exp1_pilot.csv --pi0 0.05
# confirmatory Experiment 1 (grid: configs/confirmatory_A1b.yaml once adopted in the addendum)
python scripts/04_exp1.py ... --models e2e_varnet,unet_complex,diffusion_langevin --accelerations 4,6,8 --label confirmatory --out-dir outputs/confirm
python scripts/08_stats.py --table outputs/confirm/exp1_confirmatory.csv
```
Experiments 2 to 4 reuse the same per-pair table (`04_exp1.py` records the
residuals and transfer statistics that Experiment 2 decomposes; Experiment 3
adds the shifted-acquisition conditions via `--mask-types` and the noise and
coil options in `exp3_measurement_side.py`; Experiment 4 is the phase table
built from the same rows).

## 7. Artefacts already computed (in `runs/`)

- `runs/annotated_val_A1/`: the amended bank manifest (1,725 insertions on 22
  hosts), appearance statistics, integrity manifest, and the reference
  detectability check that verified amendment A1.
- `runs/annotated_val0/`: the pre-amendment characterization (984 insertions),
  reference patches and d', the Experiment 2 and 3 demonstrations, the
  classical Experiment 1 dry run (714 pairs), and figures.
- `runs/annotated_val_A1/bank_characterization.csv` and reference patches are
  added when the CPU stage on the Mac finishes (2026-09-10 evening).

## 8. Budget (from `docs/VALIDATION.md`, grid A1b)

About 29,000 reconstructions over three models and three accelerations; the
diffusion arm is roughly 160 GPU-hours with one seed plus a 10 % replicate
subsample, VarNet and the U-Net a few hours. One rented A100-class GPU for a
week, or an institutional allocation.
