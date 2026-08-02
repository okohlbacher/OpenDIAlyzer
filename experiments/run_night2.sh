#!/bin/bash
# calib_bootstrap and calib_none both returned 0 IDs -- but they also ran at the 600 s default
# window while cirt uses its own 1435 s estimate. Two factors, one comparison. These separate them:
#
#   bootstrap + 1435 : if this recovers IDs, the WINDOW was the problem and the cheap map suffices,
#                      which would save the 424.9 s / 15277 CPU-s CiRT extraction pass outright.
#   cirt      +  600 : the reverse control. If this also collapses, window width dominates and the
#                      mapping is comparatively unimportant.
EXTRA_ARGS="-rt_calibration bootstrap -rt_extraction_window 1435" /scratch/kohlbach/bench_odia.sh calib_boot_w1435
EXTRA_ARGS="-rt_calibration cirt -rt_extraction_window 1435 -use_estimated_rt_window false" /scratch/kohlbach/bench_odia.sh calib_cirt_w1435
