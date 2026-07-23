#!/bin/sh
# Turn a DIA-NN empirical library into an OpenSWATH-ready assay library.
#
# This is the library step that nf-core pipelines route around OpenMS for:
# quantms delegates it entirely to DIA-NN, mhcquant uses EasyPQP over DDA runs.
# OpenMS itself cannot BUILD a library -- it has no RT/intensity/IM predictor,
# and OpenSwathAssayGenerator takes an existing transition list, not a FASTA.
# What OpenMS does well is REFINE one, which is what this does.
#
#   build_library.sh <diann_lib.parquet> <outdir> [swath_windows.txt]
#
# Env: OPENMS MIN_TR MAX_TR

set -e
LIB=$1
OUT=$2
WINDOWS=$3
[ -n "$OUT" ] || { echo "usage: build_library.sh <diann_lib.parquet> <outdir> [windows.txt]"; exit 2; }

OPENMS=${OPENMS:-/home/kohlbach/openms3}
export LD_LIBRARY_PATH=$OPENMS/lib:$LD_LIBRARY_PATH
HERE=$(dirname "$0")
mkdir -p "$OUT"

# OpenSWATH's default and DIA-NN's typical fragment count both land at 6.
# Keeping them equal is what makes a head-to-head comparison about the
# extraction algorithm rather than about how much evidence each tool was given.
MIN_TR=${MIN_TR:-6}
MAX_TR=${MAX_TR:-6}

echo "=== 1/3  DIA-NN parquet -> OpenSWATH TSV ==="
# Handles the two conversions that fail silently: RT minutes->seconds, and
# dropping DIA-NN's handful of decoys so OpenSWATH can generate its own.
python3 "$HERE/diann_lib_to_openswath.py" "$LIB" "$OUT/targets.tsv"

echo
echo "=== 2/3  OpenSwathAssayGenerator (refine transitions) ==="
# Selects the best MAX_TR b/y transitions per precursor and, given a windows
# file, drops fragments falling inside the precursor isolation window -- those
# are contaminated by the precursor itself and are actively misleading.
set -- -in "$OUT/targets.tsv" -in_type tsv \
       -out "$OUT/assays.tsv" -out_type tsv \
       -min_transitions "$MIN_TR" -max_transitions "$MAX_TR" \
       -allowed_fragment_types b,y -allowed_fragment_charges 1,2 \
       -threads 8
[ -n "$WINDOWS" ] && [ -f "$WINDOWS" ] && set -- "$@" -swath_windows_file "$WINDOWS"
"$OPENMS/bin/OpenSwathAssayGenerator" "$@" 2>&1 | grep -viE "^Progress|^-- done|^$" | head -10

echo
echo "=== 3/3  OpenSwathDecoyGenerator (shuffle decoys) ==="
# Decoys must come from OpenSWATH, not DIA-NN: PyProphet fits its error model
# on these, so they have to be matched to OpenSWATH's own scoring.
"$OPENMS/bin/OpenSwathDecoyGenerator" \
  -in "$OUT/assays.tsv" -in_type tsv \
  -out "$OUT/library.tsv" -out_type tsv \
  -method shuffle -threads 8 2>&1 | grep -E "Number of|Writing"

echo
echo "=== summary ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){h[$i]=i}; next}
  {n++; g[$h["TransitionGroupId"]]++; d[$h["Decoy"]]++}
  END{printf "  transitions   %d\n  precursors    %d\n  per precursor %.1f\n  targets %d  decoys %d\n",
      n, length(g), n/length(g), d[0], d[1]}' "$OUT/library.tsv"
echo "  -> $OUT/library.tsv"
