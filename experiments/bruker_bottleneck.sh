#!/usr/bin/env bash
# Is the Bruker parallelism collapse a property of the DATA or of the .d READER?
#
# OpenSwathWorkflow used ~84 cores on a TripleTOF mzML but only ~1.5 cores on a Bruker .d.
# S08 exists as BOTH a .d and an .mzpeak converted from that same .d, so holding the data,
# the library and the thread count fixed and swapping only the container isolates the reader.
#
#   .d      -> OpenMS SwathFile + opentims
#   .mzpeak -> our streaming adapter
# Same tool, same library, same threads. If %CPU stays ~1.5 cores on both, the bottleneck is
# the data/window geometry; if mzPeak scales up, it is the .d reader.
set -uo pipefail
ROOT=${ROOT:-/scratch/$USER}
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
export LD_LIBRARY_PATH=$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib:${LD_LIBRARY_PATH:-}
W=$ROOT/bench/bruker; mkdir -p "$W"
LIB=$ROOT/bench/timing/lib_1in2000.tsv          # same library both arms; identity is irrelevant here
D=$ROOT/rawd/FKL4341-S08-A-3_K30454-19-3_QKL29021A9_Slot1-12_1_1305.d
MZP=$ROOT/mzpeak_out/S08_bounds.mzpeak
NPROC=$(nproc)
log(){ printf '\n[%s] === %s ===\n' "$(date +%T)" "$*"; }

arm(){ # name input
  local name=$1 in=$2
  [[ -e $in ]] || { echo "$name: input missing ($in)"; return; }
  log "$name  <- $(basename "$in")"
  local S=$(date +%s)
  /usr/bin/time -v "$ROOT/odia-build/OpenDIAlyzer" \
    -in "$in" -tr "$LIB" -out "$W/$name.osw" \
    -threads "$NPROC" -recal_passes 1 -tempDirectory "$TMPDIR" > "$W/$name.log" 2>&1
  local RC=$? E=$(( $(date +%s) - S ))
  local CPU RSS FEAT
  CPU=$(grep -oE "Percent of CPU this job got: [0-9]+%" "$W/$name.log" | grep -oE "[0-9]+")
  RSS=$(grep -oE "Maximum resident set size \(kbytes\): [0-9]+" "$W/$name.log" | grep -oE "[0-9]+$")
  FEAT=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/$name.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
  # 0 features means the run aborted somewhere -- treat as FAILURE, never as a timing datapoint.
  local verdict=ok; [[ ${FEAT:-0} -eq 0 ]] && verdict=FAILED-0-FEATURES
  printf '    rc=%s wall=%ss cpu=%s%% (~%s cores) rss=%sMB features=%s %s\n' \
    "$RC" "$E" "${CPU:-?}" "$(( ${CPU:-0} / 100 ))" "$(( ${RSS:-0}/1024 ))" "$FEAT" "$verdict"
  echo "$name $E $RC $FEAT ${CPU:-0} ${RSS:-0}" >> "$W/results.txt"
}

rm -f "$W/results.txt"
arm bruker_d  "$D"
arm mzpeak    "$MZP"
log "results (arm wall rc features cpu_pct rss_kb)"
cat "$W/results.txt"
