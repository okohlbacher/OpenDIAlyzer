#!/usr/bin/env bash
# Stock OpenSwathWorkflow on the SAME inputs as the OpenDIAlyzer timing pilot, so the two are
# directly comparable: same .d, same library subsets, same thread count, single pass.
# This is the regression baseline -- OpenDIAlyzer must not be slower than this.
set -uo pipefail
ROOT=${ROOT:-/scratch/$USER}
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
IN=$ROOT/rawd/FKL4341-S08-A-3_K30454-19-3_QKL29021A9_Slot1-12_1_1305.d
OUT=$ROOT/bench/timing
NPROC=$(nproc)

for FRAC in 2000 500; do
  SUB=$OUT/lib_1in$FRAC.tsv
  [[ -s $SUB ]] || continue
  NPREC=$(awk -F'\t' 'NR>1{g[$23]=1} END{print length(g)}' "$SUB")
  echo "=== OSW 1-in-$FRAC: $NPREC precursors ==="
  S=$(date +%s)
  # qc:min_rsq relaxed to 0.7 to match OpenDIAlyzer: stock OpenSWATH ABORTS on a predicted
  # library at its 0.95 default, and an aborted run would report 0 features and a fake-fast time.
  /usr/bin/time -v "$ROOT/openms/bin/OpenSwathWorkflow" \
    -in "$IN" -tr "$SUB" -out_features "$OUT/osw_$FRAC.osw" \
    -threads "$NPROC" -tempDirectory "$TMPDIR" \
    -readOptions normal -Calibration:qc:min_rsq 0.7 \
    > "$OUT/osw_$FRAC.log" 2>&1
  RC=$?
  E=$(( $(date +%s) - S ))
  RSS=$(grep -oE "Maximum resident set size \(kbytes\): [0-9]+" "$OUT/osw_$FRAC.log" | grep -oE "[0-9]+$")
  FEAT=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$OUT/osw_$FRAC.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print('n/a')" 2>/dev/null)
  CPU=$(grep -oE "Percent of CPU this job got: [0-9]+%" "$OUT/osw_$FRAC.log" | grep -oE "[0-9]+%")
  echo "    rc=$RC wall=${E}s cpu=${CPU:-?} rss=$(( ${RSS:-0} / 1024 ))MB features=$FEAT"
  echo "OSW $FRAC $NPREC $E $RC $FEAT ${CPU:-?} ${RSS:-0}" >> "$OUT/baseline.txt"
  [[ $RC -ne 0 ]] && { echo "    FAILED:"; grep -iE "error|abort" "$OUT/osw_$FRAC.log" | head -4; }
done
echo "=== baseline (tool frac nprec wall rc features cpu rss_kb) ==="
cat "$OUT/baseline.txt"
