#!/bin/bash
# Stream the downloaded archives for the fastMRI+ lacunar-infarct files and measure their
# contrast as they land. Detached; safe to leave running. Usage:
#   nohup caffeinate -i -s bash scripts/chain_lacunar_extraction.sh <T9 root> [archives...] &
set -u
D="$1"; shift
ARCHIVES=("$@"); [ ${#ARCHIVES[@]} -eq 0 ] && ARCHIVES=(brain_multicoil_val_batch_1 brain_multicoil_val_batch_2 brain_multicoil_train_batch_0 brain_multicoil_val_batch_0)
HERE="$(cd "$(dirname "$0")/.." && pwd)"; export PYTHONPATH="$HERE"; cd "$HERE"
OUT="$D/extracted/lacunar"; mkdir -p "$OUT" "$D/logs"; LOG="$D/logs/lacunar_chain.log"
echo "$(date '+%F %T') start; targets: $(wc -l < "$D/logs/lacunar_infarct_files.txt") files" >> "$LOG"
for a in "${ARCHIVES[@]}"; do
  n_before=$(find "$OUT" -name "*.h5" | wc -l | tr -d ' ')
  echo "$(date '+%F %T') streaming $a" >> "$LOG"
  python3 scripts/extract_subset.py --archive "$D/archives/$a.tar.xz" --wanted "$D/logs/lacunar_infarct_files.txt" \
      --out "$OUT" --quota "AXFLAIR=28,AXT1=4" --seen-log "$D/logs/lacunar_seen_$a.csv" > "$D/logs/lacunar_extract_$a.log" 2>&1
  n_after=$(find "$OUT" -name "*.h5" | wc -l | tr -d ' ')
  echo "$(date '+%F %T') $a done: $((n_after - n_before)) new lacunar files ($n_after total)" >> "$LOG"
  if [ "$n_after" -gt 0 ]; then
    python3 scripts/02_bank.py --data "$OUT" --middle-slices 0 --stats-out "$HERE/outputs/lacunar/appearance_stats.csv" \
        --manifest-out "$HERE/outputs/lacunar/bank_manifest_unused.csv" > "$HERE/outputs/lacunar/bank_stats.log" 2>&1 || true
    grep -E "boxes|levels|focal" "$HERE/outputs/lacunar/bank_stats.log" >> "$LOG"
  fi
  [ "$n_after" -ge 32 ] && break
done
echo "$(date '+%F %T') LACUNAR CHAIN DONE" >> "$LOG"
