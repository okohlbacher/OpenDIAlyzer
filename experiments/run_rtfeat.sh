#!/bin/bash
# RT position + sqrt-deviation as scoring features. Targets the 1,261 scoring-side misses.
#
# ACCEPTANCE IS NOT THE ID COUNT ALONE. Gradient position is label-independent only if decoy
# predicted RTs mirror target ones; if they do not, the feature is a label proxy and inflates
# apparent separation without finding real peptides. DIA-NN withholds its position feature from the
# LINEAR classifier (diann.cpp:6644), which suggests they hit exactly this.
# So: compare IDs AND the target/decoy score distributions against the control.
EXTRA_ARGS="-rt_features true" /scratch/kohlbach/bench_odia.sh rtfeat_on
