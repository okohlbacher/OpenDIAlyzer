#!/usr/bin/env bash
# THREE-WAY: DIA-NN vs OpenSWATH vs OpenDIAlyzer on one DIA run, one shared library.
#
# Fairness rules:
#  * ALL THREE get the SAME library. DIA-NN is not allowed to run library-free (predict from
#    FASTA), because OSW/ODIA cannot, and that would compare search strategies rather than tools.
#  * Same thread count, same machine, run SEQUENTIALLY -- concurrent arms contaminated an
#    earlier measurement badly (OpenSWATH read 84 cores under load vs 156 idle).
#  * features == 0 is a FAILURE, never a timing datapoint (three silent aborts so far).
set -uo pipefail
ROOT=/scratch/$USER
source "$ROOT/odia-env.sh"
export HOME=$ROOT/home TMPDIR=$ROOT/tmp
export LD_LIBRARY_PATH=$ROOT/mzpeak-cpp/buildfix:$ROOT/mzpenv/lib:${LD_LIBRARY_PATH:-}
B=/ceph/ibmi/abi/oliver
IN=${IN:?set IN to the mzML}
LIB=${LIB:?set LIB}
W=$ROOT/bench/threeway; mkdir -p "$W"
N=$(nproc)
log(){ printf '\n[%s] === %s ===\n' "$(date +%T)" "$*"; }
rec(){ printf '%-14s wall=%-7s rc=%-4s ids=%-10s cpu=%-8s rss=%s MB\n' "$@" | tee -a "$W/results.txt"; }

log "input $(du -h "$IN"|cut -f1)  library $(wc -l < "$LIB") transitions  threads $N"

# ---- 1. DIA-NN ---------------------------------------------------------------
log "DIA-NN"
S=$(date +%s)
LD_LIBRARY_PATH=$B/opt/diann/diann-2.0:$LD_LIBRARY_PATH /usr/bin/time -v \
  "$B/bin/diann" --f "$IN" --lib "$LIB" --threads "$N" --qvalue 0.01 \
  --out "$W/diann.tsv" --temp "$TMPDIR" --no-prot-inf > "$W/diann.log" 2>&1
RC=$?; E=$(( $(date +%s)-S ))
# DIA-NN 2.0 writes PARQUET by default -- an --out ending in .tsv still produces
# <stem>.parquet, so counting lines in the .tsv silently yields 0.
IDS=$(python3 -c "
import glob
f=glob.glob('$W/diann*.parquet')
if f:
    import pyarrow.parquet as pq
    t=pq.read_table(f[0]); d=t.to_pydict()
    col=next((c for c in t.schema.names if c.lower() in ('precursor.id','precursor_id')), None)
    print(len(set(d[col])) if col else t.num_rows)
else:
    import csv,sys
    try: print(sum(1 for _ in open('$W/diann.tsv'))-1)
    except Exception: print(0)
" 2>/dev/null)
rec DIA-NN "${E}s" "$RC" "${IDS:-0}" \
  "$(grep -oE 'Percent of CPU this job got: [0-9]+%' "$W/diann.log"|grep -oE '[0-9]+%')" \
  "$(( $(grep -oE 'Maximum resident set size \(kbytes\): [0-9]+' "$W/diann.log"|grep -oE '[0-9]+$'||echo 0)/1024 ))"

# ---- 2. OpenSWATH ------------------------------------------------------------
log "OpenSwathWorkflow"
S=$(date +%s)
/usr/bin/time -v "$ROOT/openms/bin/OpenSwathWorkflow" -in "$IN" -tr "$LIB" \
  -out_features "$W/osw.osw" -threads "$N" -tempDirectory "$TMPDIR" \
  -readOptions normal -force > "$W/osw.log" 2>&1
RC=$?; E=$(( $(date +%s)-S ))
F=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/osw.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
rec OpenSWATH "${E}s" "$RC" "${F:-0}" \
  "$(grep -oE 'Percent of CPU this job got: [0-9]+%' "$W/osw.log"|grep -oE '[0-9]+%')" \
  "$(( $(grep -oE 'Maximum resident set size \(kbytes\): [0-9]+' "$W/osw.log"|grep -oE '[0-9]+$'||echo 0)/1024 ))"

# ---- 3. OpenDIAlyzer ---------------------------------------------------------
log "OpenDIAlyzer"
S=$(date +%s)
/usr/bin/time -v "$ROOT/odia-build/OpenDIAlyzer" -in "$IN" -tr "$LIB" \
  -out "$W/odia.osw" -threads "$N" -recal_passes 1 -tempDirectory "$TMPDIR" > "$W/odia.log" 2>&1
RC=$?; E=$(( $(date +%s)-S ))
F=$(python3 -c "
import sqlite3
try: print(sqlite3.connect('$W/odia.osw').execute('SELECT COUNT(*) FROM FEATURE').fetchone()[0])
except Exception: print(0)" 2>/dev/null)
rec OpenDIAlyzer "${E}s" "$RC" "${F:-0}" \
  "$(grep -oE 'Percent of CPU this job got: [0-9]+%' "$W/odia.log"|grep -oE '[0-9]+%')" \
  "$(( $(grep -oE 'Maximum resident set size \(kbytes\): [0-9]+' "$W/odia.log"|grep -oE '[0-9]+$'||echo 0)/1024 ))"

log "RESULTS"; cat "$W/results.txt"
echo
echo "NOTE: DIA-NN 'ids' are precursors at 1% q-value; OSW/ODIA 'ids' are RAW FEATURES with no"
echo "FDR applied. These are NOT comparable as-is -- the OSW/ODIA outputs still need scoring."
