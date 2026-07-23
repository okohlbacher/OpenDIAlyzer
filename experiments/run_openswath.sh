#!/bin/sh
# Run a properly parameterised OpenSWATH pipeline on immunopeptidomics DIA.
#
# OpenSWATH's defaults are tuned for tryptic SRM-style work and are wrong here in
# ways that would quietly lose it the comparison. Every non-default below is
# justified; an unparameterised run is not a fair baseline.
#
#   run_openswath.sh <run.mzML> <library_osw.tsv> <outdir>

set -e
MZML=$1
LIB=$2
OUT=$3
[ -n "$OUT" ] || { echo "usage: run_openswath.sh <run.mzML> <lib.tsv> <outdir>"; exit 2; }

OPENMS=${OPENMS:-/home/kohlbach/openms3}
export LD_LIBRARY_PATH=$OPENMS/lib:$LD_LIBRARY_PATH
mkdir -p "$OUT"

# -threads: default is 1. Reporting single-threaded OpenSWATH against
#  128-thread DIA-NN would be a meaningless speed comparison. Not set higher
#  because mzML ingest is serial and oversubscription measurably *hurts*
#  (224 threads was 2x slower than 1 on this machine).
THREADS=${THREADS:-32}

# -mz_extraction_window: default 50 ppm is far too loose for Orbitrap and
#  admits interference. 25 ppm is what DIA-NN's own auto-calibration measured on
#  this run (22 MS1 / 25 MS2), so both tools see the same effective tolerance.
MZ_WIN=${MZ_WIN:-25}

# -rt_extraction_window: the single most important parameter for specificity.
#  Measured on this data: candidate space at +/-300 s is ~350x less constrained
#  than at +/-10 s. The library here is EMPIRICAL (derived from a search of this
#  same run), so its RT is already in the run's own scale and a tight window is
#  legitimate -- no iRT calibrants needed. Widen for a predicted library.
RT_WIN=${RT_WIN:-120}

# This dataset has no ion mobility (confirmed by odia-info), so -1 is correct
# rather than merely a default.
$OPENMS/bin/OpenSwathWorkflow \
  -in "$MZML" \
  -tr "$LIB" \
  -out_features "$OUT/features.osw" \
  -threads "$THREADS" \
  -mz_extraction_window "$MZ_WIN" \
  -mz_extraction_window_unit ppm \
  -mz_extraction_window_ms1 "$MZ_WIN" \
  -mz_extraction_window_ms1_unit ppm \
  -rt_extraction_window "$RT_WIN" \
  -ion_mobility_window -1 \
  -tempDirectory "$OUT/tmp" \
  -readOptions cacheWorkingInMemory \
  2>&1 | tee "$OUT/openswath.log"

echo "features: $OUT/features.osw"
