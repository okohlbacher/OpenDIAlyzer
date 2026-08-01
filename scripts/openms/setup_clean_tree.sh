#!/bin/bash
# Establish a clean ODIA source tree on the compute node and lock OpenMS.
#
# Why this exists: `ext/` is gitignored, so an edit to the OpenMS core leaves no trace in this
# repository. That is how a set of ad-hoc OpenMS modifications accumulated undetected. The rules
# this script enforces:
#
#   1. vendored-patches/OpenMS/opendialyzer-openswath.patch IS the OpenMS delta. Regenerate it with
#      `git -C ext/OpenMS diff -- . ':(exclude)*PEPTDEEP*' > \
#         vendored-patches/OpenMS/opendialyzer-openswath.patch`
#      Never edit ext/ and leave it undocumented. (PEPTDEEP is covered by peptdeep-mod-support.patch.)
#   2. The installed OpenMS and its source tree are chmod a-w after a build. Modifying the core
#      then requires `chmod -R u+w`, which is a deliberate act.
#   3. ODIA's own source is the ONLY tree that is writable. All ODIA code lives here.
set -euo pipefail
NODE=${NODE:-data}
JUMP=${JUMP:-sshgw}
ROOT=${ROOT:-/scratch/kohlbach}
LOCAL=${LOCAL:-$(cd "$(dirname "$0")/../.." && pwd)}

echo "== syncing ODIA source (this tree only) =="
rsync -a --delete -e "ssh -o BatchMode=yes -J $JUMP" \
  --exclude 'ext/' --exclude 'build/' --exclude '.git/' \
  "$LOCAL/" "$NODE:$ROOT/src/OpenDIAlyzer/"

echo "== verifying OpenMS is locked =="
ssh -o BatchMode=yes -J "$JUMP" "$NODE" "
  for d in $ROOT/openms $ROOT/src/OpenMS/src; do
    if [ -w \"\$d\" ]; then echo \"WARNING: \$d is WRITABLE -- run chmod -R a-w \$d\"; else echo \"locked: \$d\"; fi
  done"

echo "== building ODIA =="
ssh -o BatchMode=yes -J "$JUMP" "$NODE" "
  export HOME=$ROOT/home TMPDIR=$ROOT/tmp
  export LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:\${LD_LIBRARY_PATH:-}
  cd $ROOT/odia-build
  # A trailing pipe swallows the build exit status, so a FAILED build used to fall through to a
  # selftest of the PREVIOUS binary and print 'selftest OK'. That is how a broken build gets
  # benchmarked. pipefail + an explicit check stops that.
  set -o pipefail
  if ! cmake --build . --target OpenDIAlyzer -j 16 2>&1 | tail -5; then
    echo 'BUILD FAILED -- not running the selftest (it would test the previous binary)' >&2
    exit 1
  fi
  cd $ROOT && ./odia-build/OpenDIAlyzer -selftest 2>&1 | grep -iE 'selftest (OK|FAILED)'"
