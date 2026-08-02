#!/bin/bash
# The tuned settings are now DEFAULTS. This must reproduce recal600's 6930 IDs with no flags:
# the derived pass-1 window comes from CiRT's own estimate (1435 s, previously rejected by the
# yield gate and the narrow-only rule), and pass 2 defaults to 600 s.
/scratch/kohlbach/bench_odia.sh tuned_default
