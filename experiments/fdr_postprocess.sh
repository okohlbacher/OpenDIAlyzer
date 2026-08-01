#!/usr/bin/env bash
# Turn raw OSW/ODIA features into IDs at 1% FDR, so they are comparable with DIA-NN's
# q-value-filtered precursor count. Without this step the three-way is meaningless:
# DIA-NN reports precursors at 1% q, the other two report every raw feature.
#
# Both .osw files get the SAME treatment (same pyprophet, same FDR level, same context),
# so any difference is the tool, not the statistics.
set -uo pipefail
ROOT=/scratch/$USER
B=/ceph/ibmi/abi/oliver
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
export OPENMS_DATA_PATH=$ROOT/openms/share/OpenMS      # pyopenms warns without it
W=$ROOT/bench/threeway
PP=$B/bin/pyprophet
log(){ printf '[%s] %s\n' "$(date +%T)" "$*"; }

for tool in osw odia; do
  f=$W/$tool.osw
  [[ -s $f ]] || { log "skip $tool (no .osw)"; continue; }
  cp "$f" "$W/$tool.scored.osw"
  log "$tool: pyprophet score"
  # --level=ms2 scores the MS2 feature level; classifier LDA matches what OpenDIAlyzer does
  # in-process, so the two are on the same statistical footing.
  "$PP" score --in "$W/$tool.scored.osw" --level=ms2 --classifier=LDA \
    > "$W/$tool.pyprophet.log" 2>&1 || { log "$tool: score FAILED"; tail -5 "$W/$tool.pyprophet.log"; continue; }
  n=$(python3 -c "
import sqlite3
c=sqlite3.connect('$W/$tool.scored.osw')
try:
    print(c.execute('''SELECT COUNT(DISTINCT f.PRECURSOR_ID) FROM FEATURE f
                       JOIN SCORE_MS2 s ON s.FEATURE_ID=f.ID
                       WHERE s.QVALUE < 0.01 AND s.RANK=1''').fetchone()[0])
except Exception as e: print('ERR',e)" 2>/dev/null)
  log "$tool: precursors at 1% FDR = $n"
  echo "$tool 1pct_fdr_precursors=$n" >> "$W/fdr.txt"
done
echo; echo "=== FDR-filtered ==="; cat "$W/fdr.txt" 2>/dev/null
