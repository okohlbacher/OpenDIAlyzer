#!/bin/bash
# outer_loop_threads and innerBatchSize only work TOGETHER.
#
#   -innerBatchSize alone (batch2k): more batches, but the inner loop had ONE thread
#                                    (max(1, 224/-1) == 1), so it bought overhead only.
#   -outer_loop_threads alone (outer16): 14 inner threads, but ONE batch per window
#                                    (423079/150 ~ 2820 compounds vs auto batch size 10000),
#                                    so 16 outer x 1 batch = 16.4 avg cores. Worse than the
#                                    38.2 baseline.
#
# The pairing has to supply both: outer windows in flight, AND enough batches inside each window
# to occupy the inner team. ~2820 compounds per window, so batch size ~ 2820 / inner_threads.
#
#   outer 16 -> inner 14 -> need ~14 batches/window -> batch ~200
#   outer 32 -> inner  7 -> need ~7  batches/window -> batch ~400
#
# Both target 224-way. If neither moves avg cores above 38.2, the window-level parallelism is the
# real ceiling and this whole line of attack is done.
EXTRA_ARGS="-outer_loop_threads 16 -innerBatchSize 200" /scratch/kohlbach/bench_odia.sh outer16_b200
EXTRA_ARGS="-outer_loop_threads 32 -innerBatchSize 400" /scratch/kohlbach/bench_odia.sh outer32_b400

# Does the default CiRT calibration earn its cost? Measured at 15277 CPU-s for 500 anchor compounds
# -- 69% of a full pass-1 extraction's CPU for 0.12% of the compounds -- returning 70 anchor pairs
# from 3897 candidates, with the MS1 component reporting flat residuals. IDs are readable even on a
# contended node, so this is answerable tonight regardless of load.
EXTRA_ARGS="-rt_calibration bootstrap" /scratch/kohlbach/bench_odia.sh calib_bootstrap
EXTRA_ARGS="-rt_calibration none"      /scratch/kohlbach/bench_odia.sh calib_none
