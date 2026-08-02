#!/bin/bash
# Pass-2 clipping confirmed: 240 -> 400 -> 600 s gave 6552 -> 6695 -> 6930 IDs and is still
# climbing. Push until it turns over -- a wider window eventually costs more in FDR (more decoy
# candidates) than it returns in coverage, and the turnover point is the answer.
EXTRA_ARGS="-rt_extraction_window 1435 -rt_extraction_window_recal 900"  /scratch/kohlbach/bench_odia.sh recal900
EXTRA_ARGS="-rt_extraction_window 1435 -rt_extraction_window_recal 1435" /scratch/kohlbach/bench_odia.sh recal1435
