#!/usr/bin/env bash
# Build the full OpenDIAlyzer stack (deps -> OpenMS -> OpenDIAlyzer) on a bare
# IBMI node, so no single node is a prerequisite for the project.
#
#   ssh <node> 'bash -s' < experiments/setup_node.sh
#
# Everything lands in /scratch/<user> (node-local, fast) -- NEVER /home, which is
# at quota, and NEVER /ceph/ibmi/it, which is the IT department's tree. Final
# results belong in /ceph/ibmi/abi. See the `ibmi-hpc` skill.
set -euo pipefail

ROOT=${ROOT:-/scratch/$USER}
ENVDIR=$ROOT/odiaenv
SRC=$ROOT/src
NPROC=$(nproc)

# /home is at quota: redirect every tool that writes a dotdir, and set the cwd
# too -- tools that write relative to cwd escape a redirected HOME.
export HOME=$ROOT/home
export TMPDIR=$ROOT/tmp
export XDG_CACHE_HOME=$ROOT/cache
export CARGO_HOME=$HOME/.cargo RUSTUP_HOME=$HOME/.rustup
mkdir -p "$HOME" "$TMPDIR" "$XDG_CACHE_HOME" "$SRC"
cd "$ROOT"

log() { printf '\n=== %s ===\n' "$*"; }

# ---- 1. dependencies ------------------------------------------------------
MAMBA=$ROOT/micromamba/micromamba
if [[ ! -x $MAMBA ]]; then
  log "installing micromamba"
  mkdir -p "$ROOT/micromamba"
  # --strip-components=1 would drop the leading "bin/", landing the binary at the
  # root rather than in bin/; extract as-is and reference bin/micromamba.
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C "$ROOT/micromamba" bin/micromamba
  MAMBA=$ROOT/micromamba/bin/micromamba
fi

if [[ ! -d $ENVDIR ]]; then
  log "creating build env"
  # GUI is off, so Qt6 is needed only for Core.
  "$MAMBA" create -y -p "$ENVDIR" -c conda-forge \
    cmake ninja pkg-config gxx_linux-64 \
    qt6-main libboost-devel eigen xerces-c libsvm glpk coin-or-cbc highs \
    libarrow libparquet libzip hdf5 icu xz libxml2 zlib bzip2 zstd libcurl \
    sqlite onnxruntime-cpp
fi
export PATH=$ENVDIR/bin:$PATH
export CMAKE_PREFIX_PATH=$ENVDIR

# ---- 2. sources -----------------------------------------------------------
log "fetching sources"
[[ -d $SRC/OpenDIAlyzer ]] \
  || git clone --quiet https://github.com/okohlbacher/OpenDIAlyzer.git "$SRC/OpenDIAlyzer"
git -C "$SRC/OpenDIAlyzer" pull --quiet --ff-only || true

OPENMS_REF=$(sed -n 's/^\*\*Base commit:\*\* `\([0-9a-f]*\)`.*/\1/p' \
              "$SRC/OpenDIAlyzer/patches/README.md" | head -1)
: "${OPENMS_REF:?could not read base commit from patches/README.md}"

if [[ ! -d $SRC/OpenMS ]]; then
  git clone --quiet https://github.com/OpenMS/OpenMS.git "$SRC/OpenMS"
  # Standing rule: never push to OpenMS.
  git -C "$SRC/OpenMS" remote set-url --push origin DISABLED_no_push
fi
git -C "$SRC/OpenMS" checkout --quiet "$OPENMS_REF"
# Reset tracked AND untracked files: the patch adds new files, so a bare
# `checkout -- .` leaves them behind and the next apply fails on re-run.
git -C "$SRC/OpenMS" checkout --quiet -- .
git -C "$SRC/OpenMS" clean --quiet -fd
log "applying OpenDIAlyzer patches to OpenMS @ $OPENMS_REF"
git -C "$SRC/OpenMS" apply "$SRC/OpenDIAlyzer/patches/openms-opendialyzer.patch"

# ---- 3. OpenMS ------------------------------------------------------------
# Re-running is normal (the OpenDIAlyzer step below iterates), and reconfiguring
# OpenMS in place does not work: its own install prefix now holds an opentims CMake
# config that references an undefined target, and the cached Opentims_DIR points
# straight at it. Skip once installed; REBUILD_OPENMS=1 forces a rebuild.
if [[ -f $ROOT/openms/lib/cmake/OpenMS/OpenMSConfig.cmake && -z ${REBUILD_OPENMS:-} ]]; then
  log "OpenMS already installed at $ROOT/openms -- skipping (REBUILD_OPENMS=1 to force)"
else
log "building OpenMS on $NPROC cores"
cmake -S "$SRC/OpenMS" -B "$ROOT/OpenMS-build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$ROOT/openms" \
  -DCMAKE_PREFIX_PATH="$ENVDIR" \
  -DBOOST_USE_STATIC=OFF \
  `# OpenMS appends only cmake/Modules + cmake/Windows to CMAKE_MODULE_PATH, but` \
  `# FindONNXRuntime.cmake sits in cmake/ -- seed the path so WITH_ONNX resolves.` \
  -DCMAKE_MODULE_PATH="$SRC/OpenMS/cmake" \
  `# On a re-run the install prefix is already populated, and OpenMS's INSTALLED` \
  `# opentims config references an undefined target -- so a reconfigure picks up the` \
  `# broken copy instead of building it from source. Keep the prefix out of find().` \
  -DCMAKE_FIND_USE_INSTALL_PREFIX=OFF \
  -DWITH_GUI=OFF -DENABLE_TUTORIALS=OFF -DENABLE_DOCS=OFF \
  -DWITH_THERMO_RAW=OFF \
  -DENABLE_UPDATE_CHECK=OFF -DPYOPENMS=OFF \
  -DWITH_ONNX=ON \
  -DONNXRuntime_INCLUDE_DIR="$ENVDIR/include/onnxruntime/core/session" \
  -DONNXRuntime_LIBRARY="$ENVDIR/lib/libonnxruntime.so"
cmake --build "$ROOT/OpenMS-build" -j "$NPROC"
cmake --install "$ROOT/OpenMS-build"
fi

# OpenMS lists the PeptDeep headers in ML/sources.cmake with a `PEPTDEEP/` prefix on
# the PARENT directory's list, and they never reach the install tree; PeptDeepModX.h
# (added by our patch) is not listed at all. OpenDIALibGen includes them directly, so
# copy the directory across rather than reworking OpenMS's install rules.
for d in PEPTDEEP ONNX; do
  [[ -d $SRC/OpenMS/src/openms/include/OpenMS/ML/$d ]] \
    && cp -r "$SRC/OpenMS/src/openms/include/OpenMS/ML/$d" "$ROOT/openms/include/OpenMS/ML/"
done

# ---- 4. OpenDIAlyzer ------------------------------------------------------
log "building OpenDIAlyzer"
cmake -S "$SRC/OpenDIAlyzer" -B "$ROOT/odia-build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$ENVDIR" \
  -DOpenMS_DIR="$ROOT/openms/lib/cmake/OpenMS"
cmake --build "$ROOT/odia-build" -j "$NPROC"

# The binaries link against conda-provided libs (xerces-c, arrow, onnxruntime, ...)
# that are not on the default loader path. Emit an env file so every downstream run
# uses the same one instead of rediscovering this.
cat > "$ROOT/odia-env.sh" <<EOF
export LD_LIBRARY_PATH=$ENVDIR/lib:$ROOT/openms/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}
export PATH=$ROOT/odia-build:$ENVDIR/bin:\$PATH
export OPENMS_DATA_PATH=$ROOT/openms/share/OpenMS
EOF
# shellcheck source=/dev/null
source "$ROOT/odia-env.sh"

log "self-test"
"$ROOT/odia-build/OpenDIAlyzer" -selftest   # single dash: TOPP flag convention

log "DONE -- OpenDIAlyzer at $ROOT/odia-build/OpenDIAlyzer (env: $ROOT/odia-env.sh)"
