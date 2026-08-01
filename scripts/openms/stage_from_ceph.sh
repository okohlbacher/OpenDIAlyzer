#!/bin/bash
# Bring a node up from the shared Ceph copy. Node switching should be an ssh and a cp,
# not a rebuild.
#
# Rationale (2026-08-01): /scratch is NODE-LOCAL. With the toolchain built only into one node's
# /scratch, that node becomes a single point of failure -- and it failed twice in an hour (`data`
# announced a shutdown; `ibminode05` sat at 48x oversubscribed). Staging node->node also broke,
# because agent-forwarded credentials do not survive the parent ssh closing. Pulling FROM Ceph on
# the target node has no second hop and therefore no credential problem.
#
# Ceph is high-latency for many-small-op reads, so the persistent copy lives there and the
# latency-sensitive inputs are copied to /scratch before the run touches them.
set -euo pipefail
NODE=${NODE:?set NODE, e.g. NODE=spock}
JUMP=${JUMP:-sshgw}
CEPH=${CEPH:-/ceph/ibmi/abi/dont-backup/kohlbach/odia}
ROOT=${ROOT:-/scratch/kohlbach}

ssh -o BatchMode=yes -J "$JUMP" "$NODE" "
  set -euo pipefail
  mkdir -p $ROOT/{bench,tmp,home,bin}
  echo '== staging inputs Ceph -> node-local /scratch (sequential copy: Ceph does this well) =='
  rsync -a --info=stats2 $CEPH/data/ $ROOT/bench/ | tail -3
  echo '== toolchain =='
  rsync -a $CEPH/openms/ $ROOT/openms/
  rsync -a $CEPH/env/    $ROOT/odiaenv/
  rsync -a $CEPH/bin/    $ROOT/bin/
  chmod +x $ROOT/bin/* 2>/dev/null || true
  echo '== verify =='
  export LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:\${LD_LIBRARY_PATH:-}
  $ROOT/bin/OpenDIAlyzer -selftest 2>&1 | grep -iE 'selftest (OK|FAILED)'
  df -h $ROOT | tail -1
"
