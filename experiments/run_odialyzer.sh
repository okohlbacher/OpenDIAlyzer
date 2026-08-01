#!/bin/bash
# First MVP comparison run for OpenDIAlyzer (two-pass recalibration engine).
# ONE 2-pass run yields BOTH arms for a clean A/B:
#   pass-1 .osw  = single-pass OpenSWATH baseline
#   final  .osw  = recalibrated two-pass (narrow window via our trafoXML)
# Usage: run_odialyzer.sh <library.tsv> <run.d> <tag>
set -uo pipefail
LIB=$1; IN=$2; TAG=$3
O=/scratch/kohlbach/opendialyzer/bench/iter/$TAG
PY=/scratch/kohlbach/mamba/envs/odia/bin/python
ODIA=/scratch/kohlbach/opendialyzer/OpenDIAlyzer/build-s1/OpenDIAlyzer
CNT=/scratch/kohlbach/opendialyzer/OpenDIAlyzer/experiments/count_osw.py
export PATH=/scratch/kohlbach/mamba/envs/odia/bin:$PATH
rm -rf "$O"; mkdir -p "$O/tmp"

# The OpenDIAlyzer parent links the onnx OpenMS build; the child OpenSwathWorkflow
# is the openms3 build and needs openms3 libs. A wrapper fixes the child env so the
# two do not fight over LD_LIBRARY_PATH.
WRAP=$O/osw.sh
cat > "$WRAP" <<'EOF'
#!/bin/bash
export LD_LIBRARY_PATH=/home/kohlbach/openms3/lib
exec /home/kohlbach/openms3/bin/OpenSwathWorkflow "$@"
EOF
chmod +x "$WRAP"

export LD_LIBRARY_PATH=/scratch/kohlbach/opendialyzer/build/openms-onnx/lib
echo "[$(date +%T)] $TAG: OpenDIAlyzer 2-pass (lib $(du -h "$LIB"|cut -f1)) ..."
"$ODIA" -in "$IN" -tr "$LIB" -out "$O/final.osw" -recal_passes 2 -threads 224 \
  -rt_extraction_window 3600 -rt_extraction_window_recal 240 \
  -mz_extraction_window 30 -mz_extraction_window_ms1 30 -ion_mobility_window 0.047 \
  -openswath_workflow "$WRAP" -tempDirectory "$O/tmp" > "$O/odia.log" 2>&1
echo "  rc=$?"; grep -iE "pass [12]/|recalibration fit|anchors|done ->|error" "$O/odia.log" | tail -6

score_arm () {  # <name> <osw>
  local NAME=$1 OSW=$2
  if [ ! -s "$OSW" ]; then echo "$NAME: MISSING $OSW"; return; fi
  pyprophet score --in "$OSW" --level ms2 --threads 32 > "$O/pyp_$NAME.log" 2>&1
  printf "%s: " "$NAME"; "$PY" "$CNT" "$OSW"
  "$PY" - "$OSW" <<'PY'
import sqlite3,sys
c=sqlite3.connect(sys.argv[1])
q="""SELECT p.DECOY,s.SCORE FROM SCORE_MS2 s JOIN FEATURE f ON s.FEATURE_ID=f.ID
     JOIN PRECURSOR p ON f.PRECURSOR_ID=p.ID WHERE s.SCORE IS NOT NULL"""
t=[];d=[]
for dec,sc in c.execute(q): (d if dec else t).append(sc)
if t and d:
    dmax=max(d); print(f"    d-score TARGETS>decoy-max {sum(1 for x in t if x>dmax):,}")
PY
}
score_arm BASELINE "$O/final.osw.pass1.osw"
score_arm RECAL    "$O/final.osw"
echo "[$(date +%T)] ODIA_DONE $TAG"
