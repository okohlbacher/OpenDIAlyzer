#!/usr/bin/env bash
# Head-to-head on a COMPLETE, small dataset: the OpenSWATH tutorial data (PASS00779,
# M. tuberculosis, TripleTOF 5600 SWATH), with its own COMPLETE library.
#
# Why not a subset of the whole-proteome library: sampling every N-th precursor strips the
# CiRT/iRT anchor peptides, so stock OpenSWATH aborts ("Need at least 3 data points to remove
# outliers") while OpenDIAlyzer survives on its built-in anchor list. That is not a fair
# comparison -- it measures our anchor fallback, not extraction. A complete library removes
# the confound entirely, and this one ships published reference results.
set -uo pipefail
ROOT=${ROOT:-/scratch/$USER}
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
SRC=/ceph/ibmi/abi/data/dia_reference/PASS00779_openswath_tutorial
W=$ROOT/bench/pass00779
mkdir -p "$W"
NPROC=$(nproc)
log(){ printf '\n[%s] === %s ===\n' "$(date +%T)" "$*"; }

# Ceph is high-latency per op; stage the run and library to node-local scratch first.
RUN=$W/R1.mzML
LIB=$W/Mtb_library.tsv
[[ -s $RUN ]] || { log "staging run"; gunzip -c "$SRC/olgas_K121026_001_SW_Wayne_R1_d00.mzML.gz" > "$RUN"; }
[[ -s $LIB ]] || { log "staging library"; cp "$SRC/Mtb_TubercuList-R27_iRT_UPS_decoy_OpenMS21compatible.tsv" "$LIB"; }
NPREC=$(awk -F'\t' 'NR>1{g[$(NF-6)]=1} END{print length(g)}' "$LIB" 2>/dev/null)
log "run $(du -h "$RUN"|cut -f1), library $(wc -l < "$LIB") transitions"

run_arm(){ # name cmd...
  local name=$1; shift
  log "$name"
  local S=$(date +%s)
  /usr/bin/time -v "$@" > "$W/$name.log" 2>&1
  local RC=$? E=$(( $(date +%s) - S ))
  local RSS CPU FEAT
  RSS=$(grep -oE "Maximum resident set size \(kbytes\): [0-9]+" "$W/$name.log" | grep -oE "[0-9]+$")
  CPU=$(grep -oE "Percent of CPU this job got: [0-9]+%" "$W/$name.log" | grep -oE "[0-9]+%")
  FEAT=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/$name.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print('n/a')" 2>/dev/null)
  echo "    rc=$RC wall=${E}s cpu=${CPU:-?} rss=$(( ${RSS:-0}/1024 ))MB features=$FEAT"
  echo "$name $E $RC $FEAT ${CPU:-?} ${RSS:-0}" >> "$W/results.txt"
  [[ $RC -ne 0 ]] && grep -iE "error|abort" "$W/$name.log" | head -3
}

run_arm openswath "$ROOT/openms/bin/OpenSwathWorkflow" \
  -in "$RUN" -tr "$LIB" -out_features "$W/openswath.osw" \
  -threads "$NPROC" -tempDirectory "$TMPDIR" -readOptions normal

run_arm opendialyzer "$ROOT/odia-build/OpenDIAlyzer" \
  -in "$RUN" -tr "$LIB" -out "$W/opendialyzer.osw" \
  -threads "$NPROC" -recal_passes 1 -tempDirectory "$TMPDIR"

log "results (arm wall rc features cpu rss_kb)"
cat "$W/results.txt"
