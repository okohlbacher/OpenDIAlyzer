#!/bin/bash
# Two independent wins measured separately: pass-2 window 864 s (+50) and -ms1_scores (+41).
# Both now defaults except ms1_scores. Test them together -- if additive, ~7020.
/scratch/kohlbach/bench_odia.sh w864_only
EXTRA_ARGS="-ms1_scores" /scratch/kohlbach/bench_odia.sh w864_ms1
