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
  # This target builds ODIA ONLY. libOpenMS is a separate, deliberately read-only install, so an
  # edit under src/OpenMS is compiled by NOTHING here and the stale .so is published as if fresh.
  # That is not hypothetical: extraction-region instrumentation was 'published and verified' three
  # times and benchmarked for ~75 minutes while strings(libOpenMS.so) contained none of it, and the
  # absent output read as 'those regions cost nothing'. Refuse rather than ship a stale library.
  # \( ... \) is load-bearing: without it -o binds as (-name '*.cpp') OR ('*.h' AND -newer), so every
  # .cpp matches and the gate fires always. sed, not head: head closes the pipe and SIGPIPEs find,
  # which under 'set -o pipefail' exits the script 141 with no message at all.
  NEWER=\$(find $ROOT/src/OpenMS/src \\( -name '*.cpp' -o -name '*.h' \\) -newer $ROOT/openms/lib/libOpenMS.so 2>/dev/null | sed -n '1,3p')
  if [ -n \"\$NEWER\" ]; then
    echo '  OpenMS sources are NEWER than the installed libOpenMS.so:' >&2
    echo \"\$NEWER\" | sed 's/^/    /' >&2
    echo '  Rebuilding OpenMS is a deliberate act (the install is chmod a-w on purpose).' >&2
    echo '  Either revert the OpenMS-side edit, or rebuild+reinstall OpenMS explicitly.' >&2
    exit 1
  fi
  # PROVE the artifact contains what was built. Three separate silent failures on 2026-08-01 gave a
  # 'published, selftest OK' toolchain that did NOT contain the change about to be benchmarked:
  #   * cmake --install blocked by a chmod a-w lock, swallowed by >/dev/null
  #   * cp to the run node refused with 'Text file busy'
  #   * an untracked header (PeptDeepModX.h) removed by a git clean, hidden by incremental build
  # -selftest passes through all of them: it exercises RT transforms, not the code under test.
  # VERIFY_SYMBOL only proves anything if the symbol is NEW IN THIS CHANGE. Passing one that already
  # existed (e.g. ChunkedColumn::resolveRaw) makes the gate a tautology -- it reported '1 matches' on
  # a library that did not contain the change being deployed.
  if [ -n '${VERIFY_SYMBOL:-}' ]; then
    N=\$(nm -DC $ROOT/openms/lib/libOpenMS.so 2>/dev/null | grep -c '${VERIFY_SYMBOL:-__none__}' || true)
    echo \"  VERIFY_SYMBOL '${VERIFY_SYMBOL:-}' -> \$N matches\"
    if [ \"\$N\" -lt 1 ]; then echo '  ABSENT - refusing to publish' >&2; exit 1; fi
  fi
  echo '== publishing to Ceph =='
  # BOTH artifacts. Publishing only the ODIA binary is how a libOpenMS-side fix reached a
  # 'verified' publish and then failed to reach the run node: the binary was new, the library was
  # four hours old, and every check passed because they were checks on the binary.
  rsync -a --delete $ROOT/openms/ $CEPH/openms/
  cp $ROOT/odia-build/OpenDIAlyzer $CEPH/bin/OpenDIAlyzer
  ls -l --time-style=+%Y-%m-%d\ %H:%M $CEPH/bin/OpenDIAlyzer
"
echo "published. stage a node with:  NODE=<node> $(dirname "$0")/stage_from_ceph.sh"
