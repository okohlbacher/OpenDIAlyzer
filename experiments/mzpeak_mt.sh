#!/usr/bin/env bash
# Verify the per-thread mzPeak decoders: does removing the global decode mutex lift the
# adapter off 1 core, and does it still produce the same features?
set -uo pipefail
ROOT=/scratch/$USER
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
export LD_LIBRARY_PATH=$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib:${LD_LIBRARY_PATH:-}
W=$ROOT/bench/bruker; mkdir -p "$W"
S=$(date +%s)
/usr/bin/time -v "$ROOT/odia-build/OpenDIAlyzer" \
  -in "$ROOT/mzpeak_out/S08_bounds.mzpeak" \
  -tr "$ROOT/bench/timing/lib_1in2000.tsv" \
  -out "$W/mzpeak_mt.osw" -threads "$(nproc)" -recal_passes 1 -tempDirectory "$TMPDIR" \
  > "$W/mzpeak_mt.log" 2>&1
RC=$?; E=$(( $(date +%s) - S ))
CPU=$(grep -oE "Percent of CPU this job got: [0-9]+%" "$W/mzpeak_mt.log" | grep -oE "[0-9]+")
RSS=$(grep -oE "Maximum resident set size \(kbytes\): [0-9]+" "$W/mzpeak_mt.log" | grep -oE "[0-9]+$")
FEAT=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/mzpeak_mt.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
V=ok; [[ ${FEAT:-0} -eq 0 ]] && V=FAILED-0-FEATURES
printf 'mzpeak_mt rc=%s wall=%ss cpu=%s%% (~%s cores) rss=%sMB features=%s %s\n' \
  "$RC" "$E" "${CPU:-?}" "$(( ${CPU:-0}/100 ))" "$(( ${RSS:-0}/1024 ))" "$FEAT" "$V" \
  | tee -a "$W/results.txt"
