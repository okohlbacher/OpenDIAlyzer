#!/bin/sh
# Run a properly parameterised OpenSWATH pipeline on immunopeptidomics DIA,
# with settings matched to what DIA-NN actually used on the same run.
#
# OpenSWATH's defaults are tuned for tryptic SRM-style work. Left alone they
# would lose the comparison for reasons that have nothing to do with the
# algorithm: -threads defaults to 1 and -mz_extraction_window to 50 ppm. Every
# non-default below is justified, and the mass/RT tolerances are taken from
# DIA-NN's own calibration log for this run so both tools see the same windows.
#
#   run_openswath.sh <run.mzML> <library_osw.tsv> <outdir>
#
# Env overrides: OPENMS THREADS OUTER MZ_WIN MZ_WIN_MS1 RT_WIN

set -e
MZML=$1
LIB=$2
OUT=$3
[ -n "$OUT" ] || { echo "usage: run_openswath.sh <run.mzML> <lib.tsv> <outdir>"; exit 2; }

OPENMS=${OPENMS:-/home/kohlbach/openms3}
export LD_LIBRARY_PATH=$OPENMS/lib:$LD_LIBRARY_PATH
mkdir -p "$OUT"

# --- Threads -----------------------------------------------------------------
# OpenSwathWorkflow nests two OpenMP levels: an outer loop over SWATH maps and
# an inner loop over peptide batches, where the inner level gets
# total_threads / outer_loop_threads. Without -outer_loop_threads the outer loop
# alone caps parallelism at the number of SWATH windows (25 here), leaving most
# of the machine idle. Setting outer = 25 and threads = 200 gives 25 x 8.
THREADS=${THREADS:-200}
OUTER=${OUTER:-25}          # = number of isolation windows in this acquisition

# Loading is serial and OpenMP spin-waiting on a 224-core box measurably hurts
# (224 threads was 2x SLOWER than 1 for mzML ingest). PASSIVE parks idle threads
# instead of burning cores.
export OMP_WAIT_POLICY=PASSIVE
export OMP_PROC_BIND=close

# --- Tolerances, matched to DIA-NN's calibration on this run -----------------
# DIA-NN log: "Calibrating with mass accuracies 22 (MS1), 25 (MS2)".
# Use exactly those so neither tool gets a tolerance advantage. DIA-NN later
# reported "Recommended MS1 mass accuracy setting: 4.6 ppm" -- tightening BOTH
# tools to that is a legitimate follow-up, but it must be done to both.
MZ_WIN=${MZ_WIN:-25}         # MS2, ppm
MZ_WIN_MS1=${MZ_WIN_MS1:-22} # MS1, ppm

# RT window. DIA-NN reports RT in MINUTES, so its "RT window set to 17.2117" is
# +/-17.2 min for a PREDICTED library. OpenSWATH's -rt_extraction_window is a
# FULL width in SECONDS, so the equivalent is 2 * 17.2 * 60 = 2064 s.
# For an EMPIRICAL library (RT already in this run's own scale) the window
# collapses to tens of seconds -- pass RT_WIN explicitly to match whatever
# DIA-NN's log reports for that run. Getting this wrong by 60x is easy; it is
# the parameter that most affects both specificity and runtime.
RT_WIN=${RT_WIN:-2064}

echo "threads=$THREADS outer=$OUTER  mz=${MZ_WIN}ppm ms1=${MZ_WIN_MS1}ppm rt_window=${RT_WIN}s"

/usr/bin/time -v $OPENMS/bin/OpenSwathWorkflow \
  -in "$MZML" \
  -tr "$LIB" \
  -out_features "$OUT/features.osw" \
  -threads "$THREADS" \
  -outer_loop_threads "$OUTER" \
  -mz_extraction_window "$MZ_WIN" \
  -mz_extraction_window_unit ppm \
  -mz_extraction_window_ms1 "$MZ_WIN_MS1" \
  -mz_extraction_window_ms1_unit ppm \
  -rt_extraction_window "$RT_WIN" \
  -ion_mobility_window -1 \
  -tempDirectory "$OUT/tmp" \
  2>&1 | tee "$OUT/openswath.log"

echo "features: $OUT/features.osw"
grep -E "Elapsed \(wall|User time|Percent of CPU|Maximum resident" "$OUT/openswath.log" || true
