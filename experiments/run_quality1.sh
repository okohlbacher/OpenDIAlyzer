#!/bin/bash
# Two one-run quality tests, both previously blocked on things now fixed.
#
# ms1_scores: measured +68 IDs and dismissed against a +/-83 "noise floor" that was itself the
# fold-assignment bug. With deterministic scoring, +68 would be real. Never retested since.
EXTRA_ARGS="-ms1_scores" /scratch/kohlbach/bench_odia.sh ms1_on
