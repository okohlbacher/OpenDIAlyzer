#!/usr/bin/env bash
# Full-library OpenDIAlyzer run on the staged S08 diaPASEF run, plus a stock
# OpenSwathWorkflow arm for comparison.
#
#   ssh <node> 'bash -s' < experiments/run_benchmark_on_node.sh
#
# Assumes setup_node.sh and build_library_on_node.sh have run.
set -euo pipefail

ROOT=${ROOT:-/scratch/$USER}
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp     # /home is at quota

IN=${IN:-/scratch/agxt/S08_full.mzML}       # node-local: Ceph is too slow for this
LIB=${LIB:-$ROOT/bench/library/library.tsv}
OUT=${OUT:-$ROOT/bench/s08_full}
NPROC=$(nproc)
mkdir -p "$OUT"

log() { printf '\n[%s] === %s ===\n' "$(date +%T)" "$*"; }

for f in "$IN" "$LIB"; do
  [[ -f $f ]] || { echo "missing input: $f" >&2; exit 1; }
done
log "input $(du -h "$IN" | cut -f1), library $(wc -l < "$LIB") rows, $NPROC threads"

log "OpenDIAlyzer: two-pass with RT recalibration"
/usr/bin/time -v "$ROOT/odia-build/OpenDIAlyzer" \
  -in "$IN" -tr "$LIB" -out "$OUT/odia.osw" \
  -threads "$NPROC" -recal_passes 2 -tempDirectory "$TMPDIR" \
  > "$OUT/odia.log" 2>&1 || { echo "OpenDIAlyzer FAILED -- see $OUT/odia.log"; tail -30 "$OUT/odia.log"; exit 2; }
log "odia done: $(du -h "$OUT/odia.osw" 2>/dev/null | cut -f1)"

# Stock OpenSWATH on the same inputs. It aborts on predicted libraries at the
# default qc:min_rsq 0.95, so relax it to the same 0.7 OpenDIAlyzer uses --
# otherwise this arm reports zero features and the comparison is meaningless.
log "OpenSwathWorkflow baseline"
/usr/bin/time -v "$ROOT/openms/bin/OpenSwathWorkflow" \
  -in "$IN" -tr "$LIB" -out_osw "$OUT/osw.osw" \
  -threads "$NPROC" -tempDirectory "$TMPDIR" \
  -Calibration:ms1_im_calibration false \
  -RTNormalization:estimateBestPeptides false \
  -RTNormalization:NrRTBins 10 \
  -RTNormalization:MinBinsFilled 8 \
  -RTNormalization:outlierMethod none \
  -qc:min_rsq 0.7 \
  > "$OUT/osw.log" 2>&1 || { echo "OpenSwathWorkflow arm failed -- see $OUT/osw.log"; tail -30 "$OUT/osw.log"; }

log "feature counts"
for f in "$OUT/odia.osw" "$OUT/osw.osw"; do
  [[ -f $f ]] || continue
  n=$(sqlite3 "$f" "SELECT COUNT(*) FROM FEATURE;" 2>/dev/null || echo "?")
  p=$(sqlite3 "$f" "SELECT COUNT(DISTINCT PRECURSOR_ID) FROM FEATURE;" 2>/dev/null || echo "?")
  printf '  %-28s features=%s  precursors=%s\n' "$(basename "$f")" "$n" "$p"
done

log "DONE -- outputs in $OUT"
