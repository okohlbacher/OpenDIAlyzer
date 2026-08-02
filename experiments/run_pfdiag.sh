#!/bin/bash
# Diagnose which prefilter criterion removes the 923 reference IDs, and test the direction I have
# NOT tested. Relaxing loses IDs (6930 -> 5580 at 3-of-6, 6567 at top-3000 peaks), so the untested
# question is whether TIGHTENING gains them -- if the prefilter is buying IDs by suppressing noise,
# 5-of-6 should beat 4-of-6. That is the control the relaxation experiments lacked.
#
# -prefilter_out exits before extraction, so the evidence dump is ~10 min rather than 40.
export ROOT=/scratch/kohlbach
LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib \
  $ROOT/bin/OpenDIAlyzer -in $ROOT/bench/astral.mzML -tr $ROOT/bench/library_ids.oswpq \
  -threads 224 -mz_extraction_window 10 -mz_extraction_window_ms1 10 \
  -prefilter_mz_extraction_window 10 -tempDirectory $ROOT/tmp \
  -prefilter_out $ROOT/bench/pf_evidence.tsv -out $ROOT/bench/pfdiag/pfdiag.oswpq \
  > $ROOT/bench/pfdiag.log 2>&1
echo "### pfdiag exit $?"
