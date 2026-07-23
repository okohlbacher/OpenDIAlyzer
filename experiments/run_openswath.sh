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

# --- Tolerances, matched to DIA-NN's OPTIMISED values on this run ------------
# Mind the factor of 2. OpenSWATH's -mz_extraction_window is a FULL width:
#     left  = mz - mz * mz_extraction_window / 2.0 * 1e-6
#     right = mz + mz * mz_extraction_window / 2.0 * 1e-6
# (ChromatogramExtractorAlgorithm.cpp:40-41). DIA-NN's "mass accuracy N ppm" is
# a +/- tolerance. So OpenSWATH needs 2N to see the same window as DIA-NN.
#
# Use the values DIA-NN OPTIMISED to, not the ones it started calibration with:
#   "[69:39] Optimised mass accuracy: 7 ppm"            -> MS2 +/-7  -> 14
#   "[10:43] Recommended MS1 mass accuracy setting: 4.6" -> MS1 +/-4.6 -> 9.2
# Its initial 22/25 were just the starting point for optimisation; matching
# those would hand OpenSWATH a ~3.5x looser tolerance and inflate its false
# positives, which would flatter neither tool honestly.
MZ_WIN=${MZ_WIN:-14}          # MS2, ppm FULL width  = 2 x DIA-NN's +/-7
MZ_WIN_MS1=${MZ_WIN_MS1:-9.2} # MS1, ppm FULL width  = 2 x DIA-NN's +/-4.6

# RT window, same full-width convention ("a value of 600 means +/- 300 s").
# DIA-NN reports RT in MINUTES, so its "RT window set to 17.2117" is +/-17.2 min
# for a PREDICTED library => full width 2 * 17.2 * 60 = 2064 s.
# For an EMPIRICAL library (RT already in this run's own scale) the window
# collapses to tens of seconds -- pass RT_WIN explicitly to match whatever
# DIA-NN's log reports for that run. Getting this wrong by 60x is easy; it is
# the parameter that most affects both specificity and runtime.
RT_WIN=${RT_WIN:-2064}

# Ion mobility. -1 means "no IM / extract over the whole range" and is correct
# for non-IM data such as PXD034539. For diaPASEF it must be set, or the IM
# dimension -- the main reason diaPASEF separates interference at all -- is
# thrown away. Use the width DIA-NN measured on the same runs (it logs
# "IM window set to ..." in 1/K0); 0.047 for the agxt S08/S23/S30 set.
IM_WIN=${IM_WIN:--1}

# -force below is required because this acquisition's isolation windows abut
# with a ~0.011 Th gap (375.43-399.43, then 399.44-423.44 -- visible in
# odia-info output). OpenSWATH aborts on ANY gap in the extraction windows.
# Here it is a rounding artefact of the method definition, ~0.05% of a 24 Th
# window, and ignoring it is correct. Do NOT keep -force for a run with
# genuinely discontinuous windows: it would silently skip real m/z ranges.

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
  -ion_mobility_window "$IM_WIN" \
  -force \
  -tempDirectory "$OUT/tmp" \
  2>&1 | tee "$OUT/openswath.log"

echo "features: $OUT/features.osw"
grep -E "Elapsed \(wall|User time|Percent of CPU|Maximum resident" "$OUT/openswath.log" || true

# OpenSwathWorkflow emits features with scores but no q-values. Without this
# step there is nothing to compare against DIA-NN's 1% FDR output -- raw feature
# counts are not identifications.
if command -v pyprophet >/dev/null 2>&1; then
  echo "--- pyprophet ---"
  pyprophet score --in "$OUT/features.osw" --level ms2 2>&1 | tail -5
else
  echo "WARNING: pyprophet not on PATH; $OUT/features.osw has no q-values yet."
  echo "         Run: micromamba run -n odia pyprophet score --in $OUT/features.osw --level ms2"
fi
