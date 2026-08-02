#!/bin/bash
# CompactLibrary on the production load path, default ON. Expected: library_load RSS drops well
# below +32.65 GB because the ~471M per-row string allocations become synthetic SSO ids, and the
# whole profile shifts down since the library is still resident when the parquet bundle is written.
# Control: same build with -compact_library false.
/scratch/kohlbach/bench_odia.sh compact_on
EXTRA_ARGS="-compact_library false" /scratch/kohlbach/bench_odia.sh compact_off
