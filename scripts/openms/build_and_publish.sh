#!/bin/bash
# Build ODIA where the source lives, then publish the binary to Ceph so any node can pick it up.
#
# The split exists because the build node and the measurement node are not the same machine and
# should not have to be: `data` carries the source tree and OpenMS build, while measurements need a
# node that is actually idle (see scripts/ibmi-nodes.sh -- on 2026-08-01 `data` sat at load 184
# while `spock` was at 0.00). Publishing through Ceph means switching is a stage, not a rebuild.
set -euo pipefail
BUILD_NODE=${BUILD_NODE:-data}
JUMP=${JUMP:-sshgw}
ROOT=${ROOT:-/scratch/kohlbach}
CEPH=${CEPH:-/ceph/ibmi/abi/dont-backup/kohlbach/odia}
LOCAL=${LOCAL:-$(cd "$(dirname "$0")/../.." && pwd)}

echo "== syncing ODIA source to $BUILD_NODE =="
rsync -a --delete -e "ssh -o BatchMode=yes -J $JUMP" \
  --exclude 'ext/' --exclude 'build/' --exclude '.git/' \
  "$LOCAL/" "$BUILD_NODE:$ROOT/src/OpenDIAlyzer/"

ssh -o BatchMode=yes -J "$JUMP" "$BUILD_NODE" "
  set -euo pipefail
  export HOME=$ROOT/home TMPDIR=$ROOT/tmp
  export LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:\${LD_LIBRARY_PATH:-}
  cd $ROOT/odia-build
  set -o pipefail
  # A trailing pipe swallows the build status; without this an failed build falls through and the
  # PREVIOUS binary gets published and benchmarked.
  if ! cmake --build . --target OpenDIAlyzer -j 16 2>&1 | tail -4; then
    echo 'BUILD FAILED -- nothing published' >&2; exit 1
  fi
  cd $ROOT && ./odia-build/OpenDIAlyzer -selftest 2>&1 | grep -iE 'selftest (OK|FAILED)'
  echo '== publishing to Ceph =='
  cp $ROOT/odia-build/OpenDIAlyzer $CEPH/bin/OpenDIAlyzer
  ls -l --time-style=+%Y-%m-%d\ %H:%M $CEPH/bin/OpenDIAlyzer
"
echo "published. stage a node with:  NODE=<node> $(dirname "$0")/stage_from_ceph.sh"
