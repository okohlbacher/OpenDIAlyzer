#!/bin/bash
# THE ASYMMETRY TEST. The prefilter scores TARGETS only; decoys enter the search because their
# TARGET passed, without their own fragment evidence ever being examined. If hit count correlates
# with score -- and it must -- searched targets carry a structural advantage unrelated to being
# real, which is anti-conservative.
#
# In the entrapment library 177,763 decoys are relabelled as targets, so they DO get evidence-tested.
# Real targets survive at 5.9% (212,292/3,603,425). If entrapment targets survive at a much lower
# rate, decoys entering the ordinary search unfiltered are systematically weaker than the targets
# they are meant to model, and the null is not comparable.
export ROOT=/scratch/kohlbach
LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib \
  $ROOT/bin/OpenDIAlyzer -in $ROOT/bench/astral.mzML -tr $ROOT/bench/library_entrap.oswpq \
  -threads 224 -mz_extraction_window 10 -mz_extraction_window_ms1 10 \
  -prefilter_mz_extraction_window 10 -tempDirectory $ROOT/tmp \
  -prefilter_out $ROOT/bench/pf_evidence_entrap.tsv -out $ROOT/bench/pfdiag/entrap_pf.oswpq \
  > $ROOT/bench/entrap_pf.log 2>&1
echo "### entrap_pf exit $?"
