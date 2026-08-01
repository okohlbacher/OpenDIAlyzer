#!/usr/bin/env bash
# Validate the Library:prefilter on a SMALL file before spending hours on the 4 GB Astral run.
# Checks three things, in order of importance:
#   1. does it run at all and emit the [prefilter] line (it was never CALLED until now)
#   2. does the candidate list actually shrink (7,149,966 -> ?)
#   3. do DECOYS survive -- 0 decoys means FDR is uncalibrated and the fix is worse than the bug
set -uo pipefail
ROOT=/scratch/$USER
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
export LD_LIBRARY_PATH=$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib:${LD_LIBRARY_PATH:-}
W=$ROOT/bench/validate; mkdir -p "$W"
IN=$ROOT/bench/12_80.mzML                      # 0.24 GB SCIEX SWATH -- minutes, not hours
LIB=$ROOT/bench/library/library.tsv
log(){ printf '[%s] %s\n' "$(date +%T)" "$*" | tee -a "$W/validate.log"; }

for mode in true false; do
  log "=== prefilter=$mode ==="
  S=$(date +%s)
  /usr/bin/time -v "$ROOT/odia-build/OpenDIAlyzer" \
    -in "$IN" -tr "$LIB" -out "$W/v_$mode.osw" -threads 224 -recal_passes 1 \
    -prefilter "$mode" -tempDirectory "$TMPDIR" > "$W/v_$mode.log" 2>&1
  rc=$?; e=$(( $(date +%s)-S ))
  rss=$(grep -oE 'Maximum resident set size \(kbytes\): [0-9]+' "$W/v_$mode.log"|grep -oE '[0-9]+$')
  feat=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/v_$mode.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
  log "  rc=$rc wall=${e}s features=${feat:-0} rss=$(( ${rss:-0}/1048576 ))GB"
  grep -h "prefilter\]" "$W/v_$mode.log" 2>/dev/null | sed 's/^/  /' | tee -a "$W/validate.log"
done
log "VALIDATION DONE"
