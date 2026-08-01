#!/bin/bash
# Resume the fully-open arm from an already-built library: OpenSWATH x3 ->
# pyprophet -> counts. S08 runs first as a fail-fast probe.
# Usage: run_searches.sh <library.tsv>
set -uo pipefail
LIB=$1
F=/scratch/kohlbach/opendialyzer/bench/agxt_full
RAW=/scratch/agxt/raw
RUNSH=/scratch/kohlbach/opendialyzer/OpenDIAlyzer/experiments/run_openswath.sh
CNT=/scratch/kohlbach/opendialyzer/OpenDIAlyzer/experiments/count_osw.py
PY=/scratch/kohlbach/mamba/envs/odia/bin/python
export PATH=/scratch/kohlbach/mamba/envs/odia/bin:$PATH
mkdir -p $F/osw

declare -A RAWF=(
  [S08]=FKL4341-S08-A-3_K30454-19-3_QKL29021A9_Slot1-12_1_1305.d
  [S23]=FKL4341-S23-A-8_K253-09-2_QKL29026AF_Slot1-17_1_1320.d
  [S30]=FKL4341-S30-A-10_K25423-16-1_QKL29028AV_Slot1-19_1_1326.d )

# The first S08 attempt aborted with "rsq: 0.866 is below limit of 0.95" --
# OpenSWATH's RT-normalisation quality gate, not a data problem (99/113 CiRT
# anchors found, 2964/2964 chromatograms non-empty). A PREDICTED library cannot
# reach r^2 0.95 against observed RT. So relax the gate, and switch alignment to
# lowess, which the docs recommend for exactly this case ("for many, noisy anchor
# points"); estimateBestPeptides drops poorly-shaped anchors before the fit.
CAL="-Calibration:qc:min_rsq 0.7 -Calibration:qc:min_coverage 0.3"
CAL="$CAL -Calibration:RTNormalization:alignmentMethod lowess"
CAL="$CAL -Calibration:RTNormalization:estimateBestPeptides"

# SCHEDULER. Setting -outer_loop_threads >= 0 requests the LEGACY nested-OpenMP
# per-SWATH path and *disables* the OpenMS 3.6 SWATH wave scheduler. That is what
# made the first attempts pathological: a predicted library is wildly imbalanced
# across SWATH windows (SWATH 2 held 1.54M of 7.2M compounds, 154 batches, vs 14
# for others), the legacy path cannot rebalance across the outer loop, and the
# run became hostage to the fattest window (~13 h critical path).
#
# The wave scheduler needs all three: outer_loop_threads -1, batchSize 0, and
# in-memory reads. Measured A/B on the mini library, same file, same calibration:
#   legacy (outer=4, threads=32) : 20:27
#   wave   (auto, 224 threads)   :  7:36   <- 2.7x, and it fixes the imbalance
SCHED="-batchSize 0 -readOptions cacheWorkingInMemory"

run_one () {
  local t=$1
  rm -rf "$F/osw/$t"
  EXTRA="$CAL $SCHED" IM_WIN=0.047 MZ_WIN=30 MZ_WIN_MS1=30 RT_WIN=3600 THREADS=224 OUTER=-1 \
    "$RUNSH" "$RAW/${RAWF[$t]}" "$LIB" "$F/osw/$t" > "$F/osw/$t.runlog" 2>&1
}

echo "[$(date +%T)] library: $LIB ($(du -h "$LIB" | cut -f1))"
echo "[$(date +%T)] S08 probe ..."
run_one S08
if ! $PY -c "import sqlite3,sys; sqlite3.connect(sys.argv[1]).execute('SELECT 1 FROM SCORE_MS2 LIMIT 1')" \
     "$F/osw/S08/features.osw" 2>/dev/null; then
  echo "[$(date +%T)] S08 has no SCORE_MS2 -- extraction or pyprophet failed:"
  grep -iE "error|rsq|exception" "$F/osw/S08/openswath.log" | head -5
  exit 2
fi
$PY "$CNT" "$F/osw/S08/features.osw" | tee "$F/osw/S08.count"
N=$($PY "$CNT" "$F/osw/S08/features.osw" | grep -oE "prec=[0-9]+" | cut -d= -f2)
echo "[$(date +%T)] S08 precursors at 1% = ${N:-0}"
[ "${N:-0}" -lt 3000 ] && { echo "[$(date +%T)] S08 LOW -- stopping."; exit 3; }

echo "[$(date +%T)] S23 then S30 (sequential, full machine each) ..."
run_one S23
run_one S30
$PY "$CNT" "$F/osw/S23/features.osw" "$F/osw/S30/features.osw" | tee "$F/osw/counts.txt"
echo "[$(date +%T)] ARM COMPLETE"
