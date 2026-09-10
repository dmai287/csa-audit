#!/bin/bash
# Run every CPU-only stage in order on an extracted fastMRI directory.
# Resumable: each stage skips work already present in its output.
# Usage: bash scripts/run_cpu_stage.sh <extracted_dir> <run_name> [workers]
# Logs: outputs/<run_name>/logs/*.log ; outputs: outputs/<run_name>/*.csv
set -u
DATA="$1"; RUN="$2"; WORKERS="${3:-4}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$HERE/outputs/$RUN"; mkdir -p "$OUT/logs"
export PYTHONPATH="$HERE"
cd "$HERE"
log(){ echo "$(date '+%F %T') $*" | tee -a "$OUT/logs/stage.log"; }
run(){ name="$1"; shift; log "=== $name: start"; python3 "$@" > "$OUT/logs/$name.log" 2>&1; rc=$?; log "=== $name: exit $rc"; return $rc; }

run 01_integrity scripts/01_integrity.py --data "$DATA" --out "$OUT/integrity.csv" || exit 1
run 02_bank scripts/02_bank.py --data "$DATA" --middle-slices "${MIDDLE_SLICES:-6}" \
    --default-contrasts "${DEFAULT_CONTRASTS:-0.3,0.6,1.0}" ${CONTRAST_LEVELS:+--contrast-levels "$CONTRAST_LEVELS"} \
    --stats-out "$OUT/appearance_stats.csv" --manifest-out "$OUT/bank_manifest.csv" || exit 1
run 02b_characterize scripts/02b_characterize_bank.py --manifest "$OUT/bank_manifest.csv" --data "$DATA" \
    --out "$OUT/bank_characterization.csv" --patches-out "$OUT/reference_patches.npz" \
    --dprime-out "$OUT/reference_dprime.csv" --workers "$WORKERS" ${LIMIT_UNITS:+--limit-units $LIMIT_UNITS} || exit 1
run exp2_physics scripts/exp2_physics_demo.py --characterization "$OUT/bank_characterization.csv" --data "$DATA" \
    --out "$OUT/exp2_physics_demo.csv" --n "${EXP2_N:-24}" || exit 1
run exp3_measurement scripts/exp3_measurement_side.py --characterization "$OUT/bank_characterization.csv" --data "$DATA" \
    --out "$OUT/exp3_measurement_side.csv" --n-lesions "${EXP3_N:-12}" || exit 1
log "=== ALL CPU STAGES DONE for $RUN"
