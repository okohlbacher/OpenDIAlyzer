#!/bin/bash
# Deploy locally-edited sources to spock and verify: build + engine --selftest (exercises
# fitTrafo_/isotonic/degenerate handling against REAL OpenMS) + the in-process LDA tests
# (confirms Task A stays green on the cluster). Run this when spock maintenance ends.
# Usage: deploy_and_verify.sh   (from the repo root on the LOCAL mac)
set -uo pipefail
REMOTE=spock
RROOT=/scratch/kohlbach/opendialyzer/OpenDIAlyzer
BUILD=$RROOT/build-s1
ONNXLIB=/scratch/kohlbach/opendialyzer/build/openms-onnx/lib

echo "[$(date +%T)] rsync changed sources -> $REMOTE:$RROOT/src ..."
rsync -av -e ssh \
  src/opendialyzer.cpp src/odia_lda.h \
  src/odia_lda_test.cpp src/odia_lda_adversarial_test.cpp src/odia_lda_realdata_test.cpp \
  "$REMOTE:$RROOT/src/" || { echo "rsync FAILED (cluster still down?)"; exit 1; }

echo "[$(date +%T)] build + verify on $REMOTE ..."
ssh "$REMOTE" bash -se <<EOF
set -uo pipefail
export LD_LIBRARY_PATH=$ONNXLIB
cd $BUILD || { echo "no build dir $BUILD"; exit 1; }

echo "--- cmake build (OpenDIAlyzer + LDA tests) ---"
cmake --build . --target OpenDIAlyzer odia-lda-test odia-lda-adversarial-test -j 32 2>&1 | tail -25
BIN=$BUILD/OpenDIAlyzer
[ -x "\$BIN" ] || { echo "BUILD FAILED: no \$BIN"; exit 2; }

echo "--- engine -selftest (isotonic + degenerate fits vs real OpenMS) ---"
"\$BIN" -selftest 2>&1 | tail -6

echo "--- in-process LDA unit test ---"
./odia-lda-test 2>&1 | tail -3
echo "--- in-process LDA adversarial test (incl. P4 calibration) ---"
./odia-lda-adversarial-test 2>&1 | tail -6
echo "DEPLOY_VERIFY_DONE"
EOF
