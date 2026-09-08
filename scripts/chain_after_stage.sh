#!/bin/bash
# Wait for the CPU stage to finish, then run Experiment 1's machinery with the
# classical reconstructors (CPU dry run / "no prior" reference rows).
# Usage: nohup caffeinate -i -s bash scripts/chain_after_stage.sh <data_dir> <run_name> [workers] &
set -u
DATA="$1"; RUN="$2"; WORKERS="${3:-4}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"; OUT="$HERE/outputs/$RUN"; export PYTHONPATH="$HERE"; cd "$HERE"
echo "$(date '+%F %T') waiting for stage '$RUN' to finish" >> "$OUT/logs/chain_exp1.log"
until grep -q "ALL CPU STAGES DONE" "$OUT/logs/stage.log" 2>/dev/null; do
  if grep -q -E "exit [1-9]" "$OUT/logs/stage.log" 2>/dev/null; then echo "$(date '+%F %T') stage failed; not starting exp1" >> "$OUT/logs/chain_exp1.log"; exit 1; fi
  sleep 120
done
echo "$(date '+%F %T') stage done; launching exp1 classical dry run ($WORKERS workers)" >> "$OUT/logs/chain_exp1.log"
python3 scripts/04_exp1.py --manifest "$OUT/bank_manifest.csv" --characterization "$OUT/bank_characterization.csv" --data "$DATA" \
  --models "${EXP1_MODELS:-zero_filled,cg_sense}" --accelerations "${EXP1_R:-4,8}" --slices-per-file "${EXP1_SLICES:-2}" \
  --workers "$WORKERS" --label dryrun_classical --out-dir "$OUT" > "$OUT/logs/exp1_dryrun_classical.log" 2>&1
echo "$(date '+%F %T') exp1 dry run exit $?; EXP1 DRY RUN DONE" >> "$OUT/logs/chain_exp1.log"
