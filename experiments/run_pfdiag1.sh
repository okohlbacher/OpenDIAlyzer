#!/bin/bash
# Dump at min_fragments=1 so SUB-THRESHOLD hit counts are recorded. At the default 4 the filter
# short-circuits and every failure reports 0, which makes the losses unattributable: the measured
# distribution over 3.6M targets is 0 -> 3,482,912 then nothing at 1/2/3 then 4 -> 107,686.
# Threshold 1 records the real count for anything with any evidence at all.
export ROOT=/scratch/kohlbach
LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib \
  $ROOT/bin/OpenDIAlyzer -in $ROOT/bench/astral.mzML -tr $ROOT/bench/library_ids.oswpq \
  -threads 224 -mz_extraction_window 10 -mz_extraction_window_ms1 10 \
  -prefilter_mz_extraction_window 10 -prefilter_min_fragments 1 -tempDirectory $ROOT/tmp \
  -prefilter_out $ROOT/bench/pf_evidence_t1.tsv -out $ROOT/bench/pfdiag/pfdiag1.oswpq \
  > $ROOT/bench/pfdiag1.log 2>&1
echo "### pfdiag1 exit $?"
