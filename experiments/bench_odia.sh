#!/bin/bash
# One benchmark invocation, one place. Inline heredocs launched via setsid inherit no
# environment, so two chained runs died on a missing libxerces before doing any work --
# and a run that dies at load looks exactly like a run that has not started yet.
#
# usage: bench_odia.sh <tag> [LD_PRELOAD_LIB]   (extra ODIA flags via EXTRA_ARGS)
set -euo pipefail

ROOT=${ROOT:-/scratch/kohlbach}
TAG=${1:?usage: bench_odia.sh <tag> [preload]}
PRELOAD=${2:-}
EXTRA_ARGS=${EXTRA_ARGS:-}

export LD_LIBRARY_PATH=$ROOT/odiaenv/lib:$ROOT/openms/lib:$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib:${LD_LIBRARY_PATH:-}
export OPENMS_TMPDIR=${OPENMS_TMPDIR:-$ROOT/tmp}
[ -n "$PRELOAD" ] && export LD_PRELOAD="$PRELOAD"

OUT=$ROOT/bench/$TAG
mkdir -p "$OUT" "$ROOT/tmp"
cd "$ROOT"

# Fail loudly at load time rather than 20 minutes later with an empty log.
ldd "$ROOT/bin/OpenDIAlyzer" | grep -q "not found" && { echo "FATAL: unresolved libs" >&2; ldd "$ROOT/bin/OpenDIAlyzer" | grep "not found" >&2; exit 1; }

# Record who else is on the node. Two runs of the SAME configuration came out 576 s and 824 s
# because the node acquired four other users mid-session, and that was only noticed by chance.
# A timing number without the load it was taken under is not a measurement.
{
  echo "== node conditions at start =="
  uptime
  ps -eo user,pcpu,comm --sort=-pcpu | awk 'NR<=6'
} > "$OUT/$TAG.nodeload.txt" 2>&1

# -in from /scratch (nvram), tempDirectory likewise: the zip must not unpack to a spinning fs.
/usr/bin/time -v "$ROOT/bin/OpenDIAlyzer" \
  -in bench/astral.mzML -tr bench/library_ids.oswpq \
  -threads 224 \
  -mz_extraction_window 10 -mz_extraction_window_ms1 10 -prefilter_mz_extraction_window 10 \
  -classifier gbt -tempDirectory "$ROOT/tmp" \
  $EXTRA_ARGS \
  -out "$OUT/$TAG.oswpq" > "$OUT/$TAG.oswpq.log" 2>&1
rc=$?
{ echo "== node conditions at end =="; uptime; } >> "$OUT/$TAG.nodeload.txt" 2>&1
echo "### $TAG exit $rc (extra: ${EXTRA_ARGS:-none}) load: $(uptime | sed 's/.*average: //')"
