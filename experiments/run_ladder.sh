#!/bin/bash
# The DIA-NN answers doc (diann.cpp:10446-10466) shows my "turnover at 900 s" conclusion was
# premature. DIA-NN's width search is a x1.2 UPWARD ladder that stops only after THREE consecutive
# non-improvements, capped at 10x the start. From 600 s the candidates are 720, 864, 1037 -- my
# 900 s result is ONE failure, not a stop, and 900 was not even on the ladder.
EXTRA_ARGS="-rt_extraction_window 1435 -rt_extraction_window_recal 720"  /scratch/kohlbach/bench_odia.sh ladder720
EXTRA_ARGS="-rt_extraction_window 1435 -rt_extraction_window_recal 864"  /scratch/kohlbach/bench_odia.sh ladder864
