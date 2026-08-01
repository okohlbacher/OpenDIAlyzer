#!/bin/bash
# The memory breakdown run. Ready to fire when the node returns.
#
# Produces, from ONE run:
#   [mem]           peak RSS per phase, sampled at 5 Hz, plus the global peak and its phase
#   [mem/alloc]     mallinfo2 at 5 checkpoints: in_use vs retained vs mmapped vs arena
#   [mem/component] library bytes as loaded and after prefiltering; pass_features_ size
#   [phase]         wall, cpu, avg cores and RSS in->out for all nine phases
set -euo pipefail
NODE=${NODE:-data}; JUMP=${JUMP:-sshgw}; ROOT=${ROOT:-/scratch/kohlbach}
TAG=${TAG:-membreak}

ssh -o BatchMode=yes -J "$JUMP" "$NODE" "
  export HOME=$ROOT/home
  mkdir -p $ROOT/bench/$TAG
  cd $ROOT
  BIN=$ROOT/odia-build/OpenDIAlyzer OUT=$ROOT/bench/$TAG/$TAG.oswpq \
  nohup ./bench_run_fair.sh \
    -in bench/astral.mzML -tr bench/library/library_ids.oswpq -threads 224 \
    -mz_extraction_window 10 -mz_extraction_window_ms1 10 \
    -prefilter_mz_extraction_window 10 -classifier gbt \
    -tempDirectory $ROOT/tmp > /dev/null 2>&1 &
  sleep 3; echo launched"

echo "watch with:"
echo "  ssh -J $JUMP $NODE 'grep -hE \"\\[mem|\\[phase\" $ROOT/bench/$TAG/$TAG.oswpq.log'"
