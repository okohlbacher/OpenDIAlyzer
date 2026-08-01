#!/usr/bin/env bash
# Why does OpenDIAlyzer find 0 features on TripleTOF SWATH where OpenSWATH finds 491,493?
#
# Chromatograms ARE extracted (iRT anchors 2350/2350 non-empty), so the data is reached; the
# failure is downstream. Our RT machinery is the main thing that differs from stock OSW, so
# disable it stepwise and see where features appear.
#   arm A: -rt_calibration none  + whole-RT window  -> library RT verbatim, closest to stock OSW
#   arm B: as A but wide fixed window               -> isolates window WIDTH from window CENTRE
# If A yields features, the RT map is wrong. If neither does, it is peak-picking/scoring.
set -uo pipefail
ROOT=/scratch/$USER
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
W=$ROOT/bench/pass00779
RUN=$W/R1.mzML; LIB=$W/Mtb_library.tsv
run(){ local n=$1; shift
  /usr/bin/time -v "$ROOT/odia-build/OpenDIAlyzer" -in "$RUN" -tr "$LIB" -out "$W/$n.osw" \
    -threads "$(nproc)" -recal_passes 1 -tempDirectory "$TMPDIR" "$@" > "$W/$n.log" 2>&1
  local rc=$?
  local f=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/$n.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
  echo "$n rc=$rc features=$f" | tee -a "$W/step3.txt"
}
rm -f "$W/step3.txt"
run rtnone -rt_calibration none -use_estimated_rt_window false -rt_extraction_window -1
run rtwide -rt_calibration none -use_estimated_rt_window false -rt_extraction_window 2000
