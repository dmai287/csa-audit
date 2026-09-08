#!/bin/bash
# Wait for the selective extraction to finish, then run the whole CPU stage.
# Detached from any interactive session; safe to leave overnight.
# Usage: nohup caffeinate -i -s bash scripts/chain_after_extraction.sh <extract_log> <data_dir> <run_name> [workers] &
set -u
LOG="$1"; DATA="$2"; RUN="$3"; WORKERS="${4:-4}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$HERE/outputs/$RUN"; mkdir -p "$OUT/logs"
echo "$(date '+%F %T') waiting for extraction: $LOG" >> "$OUT/logs/chain.log"
until grep -q -E "^done:|quota met" "$LOG" 2>/dev/null; do sleep 60; done
n=$(ls "$DATA"/*.h5 2>/dev/null | wc -l | tr -d ' ')
echo "$(date '+%F %T') extraction finished; $n volumes in $DATA; launching CPU stage '$RUN' with $WORKERS workers" >> "$OUT/logs/chain.log"
bash "$HERE/scripts/run_cpu_stage.sh" "$DATA" "$RUN" "$WORKERS" >> "$OUT/logs/chain.log" 2>&1
echo "$(date '+%F %T') chain complete (exit $?)" >> "$OUT/logs/chain.log"
