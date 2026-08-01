#!/bin/bash
# Fast iteration benchmark on the COMPREHENSIVE SUBSET (S08 middle-10-min RT slice,
# full library). Search -> pyprophet -> IDs@1% + the d-score TAIL separation metric
# (targets outscoring every decoy) -- the diagnostic that actually predicts IDs.
# Usage: bench_slice.sh <library.tsv> <tag>
set -uo pipefail
LIB=$1; IN=$2; TAG=$3
IN=$2
O=/scratch/kohlbach/opendialyzer/bench/iter/$TAG
PY=/scratch/kohlbach/mamba/envs/odia/bin/python
OMS=/home/kohlbach/openms3/bin
export LD_LIBRARY_PATH=/home/kohlbach/openms3/lib
export PATH=/scratch/kohlbach/mamba/envs/odia/bin:$PATH
rm -rf $O; mkdir -p $O
CAL="-Calibration:qc:min_rsq 0.7 -Calibration:qc:min_coverage 0.3 -Calibration:RTNormalization:alignmentMethod lowess -Calibration:RTNormalization:estimateBestPeptides"
SCHED="-batchSize 0 -readOptions cacheWorkingInMemory"

echo "[$(date +%T)] $TAG: search (lib $(du -h $LIB|cut -f1)) ..."
/usr/bin/time -v $OMS/OpenSwathWorkflow -in $IN -tr $LIB -out_features $O/features.osw \
  -threads 224 -outer_loop_threads -1 $SCHED \
  -mz_extraction_window 30 -mz_extraction_window_unit ppm \
  -mz_extraction_window_ms1 30 -mz_extraction_window_ms1_unit ppm \
  -rt_extraction_window 3600 -ion_mobility_window 0.047 -force \
  $CAL -tempDirectory $O/tmp > $O/osw.log 2>&1
grep -iE "Error|rsq" $O/osw.log | grep -vi warning | head -3
echo "[$(date +%T)] pyprophet ..."
pyprophet score --in $O/features.osw --level ms2 --threads 32 > $O/pyp.log 2>&1
echo "[$(date +%T)] results:"
$PY /scratch/kohlbach/opendialyzer/OpenDIAlyzer/experiments/count_osw.py $O/features.osw
$PY - "$O/features.osw" <<'PY'
import sqlite3,sys
c=sqlite3.connect(sys.argv[1])
q="""SELECT p.DECOY,s.SCORE FROM SCORE_MS2 s JOIN FEATURE f ON s.FEATURE_ID=f.ID
     JOIN PRECURSOR p ON f.PRECURSOR_ID=p.ID WHERE s.SCORE IS NOT NULL"""
t=[];d=[]
for dec,sc in c.execute(q): (d if dec else t).append(sc)
t.sort();d.sort()
if t and d:
    dmax=d[-1]; above=sum(1 for x in t if x>dmax)
    print(f"  d-score: target max {t[-1]:.2f} decoy max {dmax:.2f}  TARGETS>decoy-max {above:,}  (separation signal)")
PY
echo "[$(date +%T)] BENCH_DONE $TAG"
