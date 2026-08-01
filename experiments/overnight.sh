#!/usr/bin/env bash
# Overnight driver: works the P0/P1 priorities from docs/OpenDIAlyzer-perf-review.md.
#
# Runs everything SEQUENTIALLY and never alongside another benchmark -- concurrency produced
# two retracted conclusions today (OpenSWATH read 84 cores under load vs 156 idle).
# Every arm records wall / %CPU / peak RSS / feature count, and features==0 is a FAILURE,
# never a timing datapoint (three silent aborts today looked like fast results).
set -uo pipefail
ROOT=/scratch/$USER
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
export LD_LIBRARY_PATH=$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib:${LD_LIBRARY_PATH:-}
W=$ROOT/bench/overnight; mkdir -p "$W"
IN=$ROOT/bench/astral.mzML
TSV=$ROOT/bench/library/library.tsv
PQP=$ROOT/bench/library/library.pqp
R=$W/RESULTS.md
log(){ printf '[%s] %s\n' "$(date +%F' '%T)" "$*" | tee -a "$W/driver.log"; }

# --- wait for the in-flight three-way arms, whatever their fate ---------------
log "waiting for in-flight arms"
while pgrep -f "OpenDIAlyzer -in $IN" >/dev/null || pgrep -f "OpenSwathWorkflow -in $IN" >/dev/null; do sleep 120; done
log "machine clear"

arm(){ # arm <name> <notes> -- <cmd...>
  local name=$1 notes=$2; shift 3
  log "ARM $name"
  local S=$(date +%s)
  /usr/bin/time -v "$@" > "$W/$name.log" 2>&1
  local rc=$? e=$(( $(date +%s)-S ))
  local cpu rss feat osw="$W/$name.osw"
  cpu=$(grep -oE 'Percent of CPU this job got: [0-9]+%' "$W/$name.log"|grep -oE '[0-9]+')
  rss=$(grep -oE 'Maximum resident set size \(kbytes\): [0-9]+' "$W/$name.log"|grep -oE '[0-9]+$')
  feat=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
  local verdict=ok; [[ ${feat:-0} -eq 0 ]] && verdict='**FAILED (0 features)**'
  printf '| %s | %s s | %s | %s%% (~%s cores) | %.1f GB | %s | %s |\n' \
    "$name" "$e" "${feat:-0}" "${cpu:-?}" "$(( ${cpu:-0}/100 ))" \
    "$(echo "${rss:-0}/1048576"|bc -l)" "$rc" "$verdict" >> "$R"
  log "  $name rc=$rc wall=${e}s feat=${feat:-0} cpu=${cpu:-?}% rss=$(( ${rss:-0}/1048576 ))GB"
}

cat > "$R" <<'HDR'
# Overnight results — P0 memory / P1 library-load experiments

Thermo Astral plasma DIA (~4 GB mzML), complete human library, 224-core node, arms run
SEQUENTIALLY. Reference: DIA-NN 2.0 = 509 s, 8164 IDs @1% q, ~75 cores, 22.7 GB.
Processed using DIA-NN.

| arm | wall | features | CPU | peak RSS | rc | verdict |
|---|---|---|---|---|---|---|
HDR

# --- P2/P0: does the evidence prefilter stop the feature explosion? -----------
# Pre-fix reference on this exact input: 27,344,163 features, 1753 GB, 2.9 h, with
# OpenSwathWorkflow retaining 34,165 of 3,603,425 target precursors. If the prefilter works,
# features and RSS should both fall by ~2 orders of magnitude AND decoys must survive
# (grep the prefilter line in the log -- 0 decoys means FDR is uncalibrated and the run is void).
arm prefilter_on "Library:prefilter ON (new default)" -- "$ROOT/odia-build/OpenDIAlyzer" \
  -in "$IN" -tr "$TSV" -out "$W/prefilter_on.osw" -threads 224 -recal_passes 1 \
  -prefilter true -tempDirectory "$TMPDIR"
grep -h "prefilter\]" "$W/prefilter_on.log" 2>/dev/null | tee -a "$R" | tee -a "$W/driver.log"

# --- P1: does a PQP library remove the 7:39 serial TSV parse? -----------------
if [[ -s $PQP ]]; then
  arm p1_pqp "PQP library + prefilter" -- "$ROOT/odia-build/OpenDIAlyzer" \
    -in "$IN" -tr "$PQP" -out "$W/p1_pqp.osw" -threads 224 -recal_passes 1 \
    -prefilter true -tempDirectory "$TMPDIR"
else
  log "SKIP p1_pqp: no $PQP (conversion failed?)"
fi

# --- P0: is the 861 GB unbatched result accumulation? -------------------------
LIB=$([[ -s $PQP ]] && echo "$PQP" || echo "$TSV")
arm p0_batch5k "batchSize 5000" -- "$ROOT/odia-build/OpenDIAlyzer" \
  -in "$IN" -tr "$LIB" -out "$W/p0_batch5k.osw" -threads 224 -recal_passes 1 \
  -batchSize 5000 -tempDirectory "$TMPDIR"

arm p0_batch50k "batchSize 50000" -- "$ROOT/odia-build/OpenDIAlyzer" \
  -in "$IN" -tr "$LIB" -out "$W/p0_batch50k.osw" -threads 224 -recal_passes 1 \
  -batchSize 50000 -tempDirectory "$TMPDIR"

# --- P0: does streaming input (no in-memory SWATH map) cut RSS, and at what CPU cost?
arm p0_stream "readOptions normal" -- "$ROOT/odia-build/OpenDIAlyzer" \
  -in "$IN" -tr "$LIB" -out "$W/p0_stream.osw" -threads 224 -recal_passes 1 \
  -readOptions normal -tempDirectory "$TMPDIR"

# --- FDR: make feature counts comparable with DIA-NN's 8164 -------------------
log "pyprophet FDR on every arm that produced features"
for f in "$W"/*.osw; do
  [[ -s $f ]] || continue
  n=$(basename "$f" .osw)
  cp "$f" "$W/$n.scored.osw"
  OPENMS_DATA_PATH=$ROOT/openms/share/OpenMS /ceph/ibmi/abi/oliver/bin/pyprophet score \
    --in "$W/$n.scored.osw" --level=ms2 --classifier=LDA > "$W/$n.pyprophet.log" 2>&1 || \
    { log "  $n: pyprophet failed"; continue; }
  ids=$(python3 -c "
import sqlite3
c=sqlite3.connect('$W/$n.scored.osw')
try: print(c.execute('''SELECT COUNT(DISTINCT f.PRECURSOR_ID) FROM FEATURE f
   JOIN SCORE_MS2 s ON s.FEATURE_ID=f.ID WHERE s.QVALUE<0.01 AND s.RANK=1''').fetchone()[0])
except Exception: print('n/a')" 2>/dev/null)
  log "  $n: $ids precursors @1% FDR"
  echo "- \`$n\`: **$ids** precursors at 1% FDR (DIA-NN: 8164)" >> "$R"
done

log "DONE"
echo >> "$R"; echo "Completed $(date)." >> "$R"
