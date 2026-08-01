#!/bin/bash
# One benchmark invocation, one place. Inline heredocs launched via setsid inherit no
# environment, so two chained runs died on a missing libxerces before doing any work --
# and a run that dies at load looks exactly like a run that has not started yet.
#
# usage: bench_odia.sh <tag> [LD_PRELOAD_LIB]
set -euo pipefail

ROOT=${ROOT:-/scratch/kohlbach}
TAG=${1:?usage: bench_odia.sh <tag> [preload]}
PRELOAD=${2:-}

export LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib:${LD_LIBRARY_PATH:-}
export OPENMS_TMPDIR=${OPENMS_TMPDIR:-$ROOT/tmp}
[ -n "$PRELOAD" ] && export LD_PRELOAD="$PRELOAD"

OUT=$ROOT/bench/$TAG
mkdir -p "$OUT" "$ROOT/tmp"
cd "$ROOT"

# Fail loudly at load time rather than 20 minutes later with an empty log.
ldd "$ROOT/bin/OpenDIAlyzer" | grep -q "not found" && { echo "FATAL: unresolved libs" >&2; ldd "$ROOT/bin/OpenDIAlyzer" | grep "not found" >&2; exit 1; }

# -in from /scratch (nvram), tempDirectory likewise: the zip must not unpack to a spinning fs.
/usr/bin/time -v "$ROOT/bin/OpenDIAlyzer" \
  -in bench/astral.mzML -tr bench/library_ids.oswpq \
  -threads 224 \
  -mz_extraction_window 10 -mz_extraction_window_ms1 10 -prefilter_mz_extraction_window 10 \
  -classifier gbt -tempDirectory "$ROOT/tmp" \
  -out "$OUT/$TAG.oswpq" > "$OUT/$TAG.oswpq.log" 2>&1
echo "### $TAG exit $?"
