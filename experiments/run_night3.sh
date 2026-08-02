#!/bin/bash
# Verify the cumulative-rescale fix: bootstrap returned 0 IDs because pass 2 applied t2 to
# t1(rt_orig) while t2 was fitted on rt_orig. Pass 1 was healthy (p95 residual 112.86 s, better
# than cirt's 117.24 s), so if the fix is right this should now produce IDs.
EXTRA_ARGS="-rt_calibration bootstrap" /scratch/kohlbach/bench_odia.sh boot_fixed
EXTRA_ARGS="-rt_calibration none"      /scratch/kohlbach/bench_odia.sh none_fixed
