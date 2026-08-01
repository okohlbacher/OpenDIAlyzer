#!/bin/bash
# Sampling profiler for the extraction phase, via gdb backtraces.
#
# perf is present on these nodes but perf_event_paranoid=4 blocks unprivileged sampling, and
# extraction is entirely inside OpenMS's performExtraction() -- so adding timers means patching the
# core. gdb attaching to our own process needs no privilege and answers the question that matters:
# WHERE do extraction's 48,120 CPU-seconds go? Every proposal about extraction is a guess until
# this is known.
#
# Aggregates function names across all threads. Not exact (backtrace sampling is biased toward
# frames that are slow to unwind) but decisive for finding a dominant symbol.
set -uo pipefail
PID=${PID:?set PID}
N=${N:-25}
SLEEP=${SLEEP:-3}
OUT=${OUT:-/dev/shm/extract_profile.txt}

: > "$OUT"
for i in $(seq 1 "$N"); do
  gdb -p "$PID" -batch -ex "set pagination off" -ex "thread apply all bt 12" 2>/dev/null \
    | grep -oE '\bin [A-Za-z_][A-Za-z0-9_:<>~]*' | sed 's/^in //' >> "$OUT"
  sleep "$SLEEP"
done
echo "=== top frames across $N samples ==="
sort "$OUT" | uniq -c | sort -rn | head -25
