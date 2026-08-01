#!/usr/bin/env bash
# OpenSWATH arm with parameters MATCHED to OpenDIAlyzer, so the only variable is the tool.
# Stock OSW defaults differ: mz window 50 ppm (ODIA 30), qc:min_rsq 0.95 (ODIA 0.7 -- 0.95
# assumes spike-in iRT kits and ABORTS on a predicted library), IM window on.
set -uo pipefail
ROOT=/scratch/$USER
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
W=$ROOT/bench/threeway
S=$(date +%s)
/usr/bin/time -v "$ROOT/openms/bin/OpenSwathWorkflow" \
  -in "$ROOT/bench/astral.mzML" -tr "$ROOT/bench/library/library.tsv" \
  -out_features "$W/osw2.osw" -threads 24 -tempDirectory "$TMPDIR" \
  -readOptions normal -force \
  -Calibration:qc:min_rsq 0.7 \
  -mz_extraction_window 30 -mz_extraction_window_unit ppm \
  -mz_extraction_window_ms1 30 -mz_extraction_window_ms1_unit ppm \
  -rt_extraction_window 600 \
  -ion_mobility_window -1 \
  > "$W/osw2.log" 2>&1
RC=$?; E=$(( $(date +%s)-S ))
F=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/osw2.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
printf 'OpenSWATH      wall=%-7s rc=%-4s features=%-10s cpu=%-8s rss=%s MB\n' \
  "${E}s" "$RC" "${F:-0}" \
  "$(grep -oE 'Percent of CPU this job got: [0-9]+%' "$W/osw2.log"|grep -oE '[0-9]+%')" \
  "$(( $(grep -oE 'Maximum resident set size \(kbytes\): [0-9]+' "$W/osw2.log"|grep -oE '[0-9]+$'||echo 0)/1024 ))" \
  | tee -a "$W/results.txt"
