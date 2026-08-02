#!/bin/bash
# Two findings drive these:
#  1. The yield gate rejects CiRT's own 1435 s estimate at 1.8% anchor yield and falls back to
#     600 s, costing 122 IDs (6430 -> 6552). Passing the window explicitly bypasses it.
#  2. Pass 2 runs at rt_extraction_window_recal = 240 s (+/-120), and the final residual p99 is
#     114.8 / 113.4 s -- pressed against the window. The distribution is clipped by pass 2, which
#     is why a 2.4x wider PASS 1 changed the final residual by only 1.4 s while adding IDs: the
#     pass-1 gain showed up in the ANCHOR spread (117.2 -> 128.5 s p95), not the final one.
# If pass 2 is genuinely clipping, widening it gains IDs. If it is not, IDs fall as the extra
# candidates cost more in FDR than they return -- which is the honest counter-hypothesis.
EXTRA_ARGS="-rt_extraction_window 1435 -rt_extraction_window_recal 400" /scratch/kohlbach/bench_odia.sh recal400
EXTRA_ARGS="-rt_extraction_window 1435 -rt_extraction_window_recal 600" /scratch/kohlbach/bench_odia.sh recal600
