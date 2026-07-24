#!/bin/bash
# Fully-open arm, full agxt proteome: OpenDIALibGen(-raw) -> OpenSwathAssayGenerator
# -> OpenSwathDecoyGenerator -> OpenSWATH x3 -> pyprophet -> counts.
#
# The in-process refine_and_decoy path is bypassed: it dies with
# "UniMod:4 ... Modification not found" at proteome scale, and MRMDecoy's default
# 0.8 min_decoy_fraction is too strict for a peptidoform library. The standalone
# TOPP tools do the same work in ~1 min per 270k precursors.
# Usage: run_full_arm.sh <raw_libgen_pid>
set -uo pipefail
RAWPID=$1
F=/scratch/kohlbach/opendialyzer/bench/agxt_full
OMS=/home/kohlbach/openms3/bin
RAW=/scratch/agxt/raw
RUNSH=/scratch/kohlbach/opendialyzer/OpenDIAlyzer/experiments/run_openswath.sh
CNT=/scratch/kohlbach/opendialyzer/OpenDIAlyzer/experiments/count_osw.py
PY=/scratch/kohlbach/mamba/envs/odia/bin/python
export LD_LIBRARY_PATH=/home/kohlbach/openms3/lib
export PATH=/scratch/kohlbach/mamba/envs/odia/bin:$PATH
mkdir -p $F/osw

declare -A RAWF=(
  [S08]=FKL4341-S08-A-3_K30454-19-3_QKL29021A9_Slot1-12_1_1305.d
  [S23]=FKL4341-S23-A-8_K253-09-2_QKL29026AF_Slot1-17_1_1320.d
  [S30]=FKL4341-S30-A-10_K25423-16-1_QKL29028AV_Slot1-19_1_1326.d )

run_one () {
  local t=$1
  IM_WIN=0.047 MZ_WIN=30 MZ_WIN_MS1=30 RT_WIN=2064 THREADS=48 OUTER=25 \
    "$RUNSH" "$RAW/${RAWF[$t]}" "$F/library.tsv" "$F/osw/$t" > "$F/osw/$t.runlog" 2>&1
}

echo "[$(date +%T)] waiting for raw libgen PID $RAWPID ..."
while kill -0 "$RAWPID" 2>/dev/null; do sleep 60; done
[ -s "$F/odia_raw.tsv" ] || { echo "FATAL: no raw TSV"; exit 1; }
echo "[$(date +%T)] raw TSV: $(du -h $F/odia_raw.tsv | cut -f1)"

echo "[$(date +%T)] OpenSwathAssayGenerator ..."
$OMS/OpenSwathAssayGenerator -in "$F/odia_raw.tsv" -out "$F/assay.tsv" > "$F/assay.log" 2>&1 \
  || { echo "FATAL: assay generator (see assay.log)"; exit 2; }
echo "[$(date +%T)] OpenSwathDecoyGenerator ..."
$OMS/OpenSwathDecoyGenerator -in "$F/assay.tsv" -out "$F/library.tsv" \
  -method shuffle -switchKR true -min_decoy_fraction 0.1 > "$F/decoy.log" 2>&1 \
  || { echo "FATAL: decoy generator (see decoy.log)"; exit 3; }
echo "[$(date +%T)] library ready: $(du -h $F/library.tsv | cut -f1)"

# S08 first as a fail-fast probe: PeptDeep's 0-1 RT scale relies on OpenSWATH
# self-anchoring RT on CiRT peptides in the library. If that breaks, IDs collapse
# and there is no point burning wall-clock on S23/S30.
echo "[$(date +%T)] S08 probe ..."
run_one S08
$PY "$CNT" "$F/osw/S08/features.osw" | tee "$F/osw/S08.count"
N=$($PY "$CNT" "$F/osw/S08/features.osw" | grep -oE "prec=[0-9]+" | cut -d= -f2)
echo "[$(date +%T)] S08 precursors at 1% = ${N:-0}"
if [ "${N:-0}" -lt 3000 ]; then
  echo "[$(date +%T)] S08 LOW -> RT calibration suspect. Stopping before S23/S30."
  exit 4
fi

echo "[$(date +%T)] S23 + S30 concurrently ..."
run_one S23 & run_one S30 & wait
$PY "$CNT" "$F/osw/S23/features.osw" "$F/osw/S30/features.osw" | tee -a "$F/osw/counts.txt"
echo "[$(date +%T)] ARM COMPLETE"
