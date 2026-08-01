#!/bin/bash
# Written locally and rsynced, not heredoc'd over ssh -- backticks and $ in comments get evaluated
# by the local shell inside a double-quoted remote command, which has now silently mangled two
# scripts today.
#
# The cacheWorkingInMemory route is dropped. It reaches 83.7 cores and +69 IDs, but dia_run_load
# costs +215 s at 3.2 avg cores to build the cache, against only 139 s saved in extraction --
# disqualifying whatever the spin policy does, so isolating passive would not change the decision.
#
# These target the SAME defect on the streaming path, with no residency requirement and no caching:
# the legacy inner batch team is sized max(1, 224 / threads_outer_loop_) = 1 when ODIA passes -1,
# so the inner loop is serial and all parallelism comes from the ~150 outer SWATH windows.
EXTRA_ARGS="-outer_loop_threads 16" /scratch/kohlbach/bench_odia.sh outer16
EXTRA_ARGS="-outer_loop_threads 32" /scratch/kohlbach/bench_odia.sh outer32
EXTRA_ARGS="-outer_loop_threads 8"  /scratch/kohlbach/bench_odia.sh outer8

# Does the default CiRT calibration earn its 260.6 s (14% of wall, 1.8% anchor yield, MS1 residuals
# flagged flat by its own diagnostic)? Deterministic scoring makes this a single-run question.
EXTRA_ARGS="-rt_calibration bootstrap" /scratch/kohlbach/bench_odia.sh calib_bootstrap
EXTRA_ARGS="-rt_calibration none"      /scratch/kohlbach/bench_odia.sh calib_none
