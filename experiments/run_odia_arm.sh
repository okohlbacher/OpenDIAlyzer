#!/bin/bash
# Fully-open arm: OpenSWATH searching an OpenDIALibGen predicted library.
# Waits for the library-gen PID, runs S08 as a fail-fast probe (RT calibration
# on the PeptDeep 0-1 scale is the one real risk), then S23+S30 concurrently.
# Usage: run_odia_arm.sh <libgen_pid>
set -u
LGPID=$1
O=/scratch/kohlbach/opendialyzer/bench/agxt_arms
LIB=$O/odia_lib/library.tsv
RAW=/scratch/agxt/raw
RUNSH=/scratch/kohlbach/opendialyzer/OpenDIAlyzer/experiments/run_openswath.sh
export OPENMS=/home/kohlbach/openms3
export LD_LIBRARY_PATH=$OPENMS/lib
export PATH=/scratch/kohlbach/mamba/envs/odia/bin:$PATH
mkdir -p $O/osw_predlib

declare -A RAWF=(
  [S08]=FKL4341-S08-A-3_K30454-19-3_QKL29021A9_Slot1-12_1_1305.d
  [S23]=FKL4341-S23-A-8_K253-09-2_QKL29026AF_Slot1-17_1_1320.d
  [S30]=FKL4341-S30-A-10_K25423-16-1_QKL29028AV_Slot1-19_1_1326.d )

count_ids () { # <features.osw> -> "total target" at q<0.01
  python3 - "$1" <<'PY'
import sqlite3,sys
c=sqlite3.connect(sys.argv[1])
tot=c.execute("SELECT COUNT(*) FROM SCORE_MS2 WHERE QVALUE<0.01").fetchone()[0]
try:
    tgt=c.execute("""SELECT COUNT(*) FROM SCORE_MS2 s JOIN FEATURE f ON s.FEATURE_ID=f.ID
                     JOIN PRECURSOR p ON f.PRECURSOR_ID=p.ID
                     WHERE s.QVALUE<0.01 AND p.DECOY=0""").fetchone()[0]
except Exception: tgt=-1
print(tot,tgt)
PY
}

run_one () { # <tag>
  local t=$1
  IM_WIN=0.047 MZ_WIN=30 MZ_WIN_MS1=30 RT_WIN=2064 THREADS=48 OUTER=25 \
    "$RUNSH" "$RAW/${RAWF[$t]}" "$LIB" "$O/osw_predlib/$t" > "$O/osw_predlib/$t.runlog" 2>&1
}

echo "[$(date +%T)] waiting for libgen PID $LGPID ..."
while kill -0 "$LGPID" 2>/dev/null; do sleep 30; done
echo "[$(date +%T)] libgen finished."
[ -s "$LIB" ] || { echo "FATAL: library missing/empty: $LIB"; exit 1; }
echo "library: $(wc -l < "$LIB") lines, $(du -h "$LIB" | cut -f1)"

echo "[$(date +%T)] S08 probe ..."
run_one S08
read TOT TGT < <(count_ids "$O/osw_predlib/S08/features.osw")
echo "[$(date +%T)] S08: total=$TOT target=$TGT peakgroups at q<0.01"
if [ "${TGT:-0}" -lt 8000 ] && [ "${TOT:-0}" -lt 8000 ]; then
  echo "[$(date +%T)] S08 LOW -> RT/IM calibration likely off. Stopping before S23/S30."
  exit 3
fi

echo "[$(date +%T)] S08 sane -> S23 + S30 concurrently ..."
run_one S23 & run_one S30 & wait
for t in S23 S30; do
  read TOT TGT < <(count_ids "$O/osw_predlib/$t/features.osw")
  echo "[$(date +%T)] $t: total=$TOT target=$TGT peakgroups at q<0.01"
done
echo "[$(date +%T)] ARM COMPLETE"
