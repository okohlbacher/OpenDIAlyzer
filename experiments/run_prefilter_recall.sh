#!/bin/bash
# Prefilter recall is 88.1%: 923 of DIA-NN's 7787 IDs are deleted before extraction. The rule is
# 4-of-top-6 predicted fragments inside one spectrum's top-1000 peaks. Relax each dimension
# separately so the responsible criterion is identified rather than guessed.
EXTRA_ARGS="-prefilter_min_fragments 3" /scratch/kohlbach/bench_odia.sh pf_frag3
EXTRA_ARGS="-prefilter_top_peaks 3000"  /scratch/kohlbach/bench_odia.sh pf_peaks3k
