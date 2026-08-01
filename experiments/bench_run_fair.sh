#!/bin/bash
# Benchmark runner for a SHARED node.
#
# The node this runs on has ~18 users and a load average that has been observed between 180 and 233
# on 224 cores. Wall clock and CPU% therefore measure how many cores the run WON, not how well it
# scales. Peak RSS and total CPU-seconds are the contention-independent numbers; everything else
# needs the load context this script records on both sides of the run.
set -u
# The binary's RUNPATH resolves these only in an interactive shell; a non-login ssh command shell
# does not source the profile that sets them, and the run dies on libxerces at startup.
export LD_LIBRARY_PATH=/scratch/kohlbach/odiaenv/lib:/scratch/kohlbach/openms/lib:${LD_LIBRARY_PATH:-}
export HOME=${HOME:-/scratch/kohlbach/home}
BIN=${BIN:?set BIN}
OUT=${OUT:?set OUT}
TAG=$(basename "$OUT")
mkdir -p "$(dirname "$OUT")"
LOG="$(dirname "$OUT")/$TAG.log"

snap() {   # load context, so a slow run can be told apart from a busy node
  echo "### load $1: $(uptime | sed 's/.*load average/load/')"
  echo "### other-user CPU $1: $(ps -eo pcpu=,user= | awk '$2!="'"$USER"'"{s+=$1} END{printf "%.0f%%", s}')"
  echo "### free GB $1: $(free -g | awk '/^Mem:/{print $7}')"
}

{ snap before; } > "$LOG"
/usr/bin/time -v "$BIN" "$@" -out "$OUT" >> "$LOG" 2>&1
rc=$?
{ snap after; } >> "$LOG"
echo "### exit $rc" >> "$LOG"
exit $rc
