#!/usr/bin/env bash
# Build the predicted assay library on whatever node we are on:
#   fasta -> OpenDIALibGen (raw) -> m/z trim -> AssayGenerator -> DecoyGenerator
#
#   ssh <node> 'bash -s' < experiments/build_library_on_node.sh
#
# Assumes experiments/setup_node.sh has run (it writes $ROOT/odia-env.sh).
set -euo pipefail

ROOT=${ROOT:-/scratch/$USER}
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp     # /home is at quota

FASTA=${FASTA:-/ceph/ibmi/abi/projects/AGXT_proteogenomics/fasta/UP000005640_human_reviewed.fasta}
OUT=${OUT:-$ROOT/bench/library}
MODELS=$ROOT/openms/share/OpenMS/models
BIN=$ROOT/odia-build
OMS=$ROOT/openms/bin
NPROC=$(nproc)
mkdir -p "$OUT"

log() { printf '\n[%s] === %s ===\n' "$(date +%T)" "$*"; }

# Ceph is high-latency for many-small-op reads; stage the FASTA locally first.
LOCAL_FASTA=$OUT/$(basename "$FASTA")
[[ -f $LOCAL_FASTA ]] || cp "$FASTA" "$LOCAL_FASTA"

log "OpenDIALibGen: predicting RT/MS2 over $(basename "$FASTA") on $NPROC threads"
"$BIN/OpenDIALibGen" -fasta "$LOCAL_FASTA" -out "$OUT/raw.tsv" \
  -model_dir "$MODELS" -threads "$NPROC" \
  -max_charge 4 -max_var_mods 1 -raw \
  > "$OUT/raw_gen.log" 2>&1
log "raw: $(du -h "$OUT/raw.tsv" | cut -f1), $(wc -l < "$OUT/raw.tsv") rows"

# Match DIA-NN's search space (its library-generation setting) so precursor counts
# are comparable rather than reflecting a wider net on our side.
log "trim to precursor m/z 300-1200"
awk -F'\t' -v OFS='\t' 'NR==1{print;next}{mz=$1+0; if(mz>=300&&mz<=1200) print}' \
  "$OUT/raw.tsv" > "$OUT/raw_trim.tsv"
log "trimmed: $(wc -l < "$OUT/raw_trim.tsv") rows"

# Defaults (product floor 350, 6 transitions) gut a DIA library: floor 200 keeps
# low-mass base peaks and 6-12 transitions matches DIA-NN's count.
log "OpenSwathAssayGenerator"
"$OMS/OpenSwathAssayGenerator" -in "$OUT/raw_trim.tsv" -out "$OUT/assay.tsv" \
  -min_transitions 6 -max_transitions 12 \
  -product_lower_mz_limit 200 -product_upper_mz_limit 2000 \
  -precursor_lower_mz_limit 300 -precursor_upper_mz_limit 1200 \
  > "$OUT/assay.log" 2>&1

log "OpenSwathDecoyGenerator (shuffle)"
"$OMS/OpenSwathDecoyGenerator" -in "$OUT/assay.tsv" -out "$OUT/library.tsv" \
  -method shuffle -switchKR true -min_decoy_fraction 0.1 \
  > "$OUT/decoy.log" 2>&1

log "LIBRARY READY: $OUT/library.tsv ($(du -h "$OUT/library.tsv" | cut -f1))"
echo "target/decoy split:"
awk -F'\t' 'NR>1{print $25}' "$OUT/library.tsv" | sort | uniq -c
