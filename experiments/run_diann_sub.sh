#!/bin/bash
# DIA-NN arm on the SAME protein subset OpenDIALibGen used, so both tools face an
# identical search space. Digest/mod params copied from the in-house predicted-library
# gen log (cuts K*,R*, MC1, len 7-30, m/z 300-1200, z2-4, Cam fixed, 1 var Ox(M));
# search tolerances copied from run-odia.sh so this matches the existing arms.
set -uo pipefail
FASTA=${1:-/scratch/kohlbach/opendialyzer/tmp/sub1000.fasta}
OUT=${2:-/scratch/kohlbach/opendialyzer/bench/inhouse_sub1000/diann}
D=/home/kohlbach/diann
mkdir -p "$OUT"
ARGS=()
for f in /scratch/inhouse/raw/IH1.d /scratch/inhouse/raw/IH3.d /scratch/inhouse/raw/IH2.d; do
  ARGS+=(--f "$f")
done
cd "$D"
export LD_LIBRARY_PATH=$D
echo "==== START $(date) ===="
/usr/bin/time -v ./diann-linux \
  "${ARGS[@]}" \
  --fasta "$FASTA" \
  --fasta-search --gen-spec-lib --predictor \
  --out-lib "$OUT/lib.predicted.speclib" \
  --out "$OUT/report.parquet" \
  --threads 96 \
  --qvalue 0.5 \
  --cut "K*,R*" \
  --missed-cleavages 1 \
  --min-pep-len 7 --max-pep-len 30 \
  --min-pr-mz 300 --max-pr-mz 1200 \
  --min-pr-charge 2 --max-pr-charge 4 \
  --unimod4 --var-mods 1 --var-mod "UniMod:35,15.9949,M" \
  --mass-acc 7 --mass-acc-ms1 10 --window 6 \
  --verbose 1 > "$OUT/run.log" 2>&1
echo "EXIT=$?"
echo "==== DONE $(date) ===="
