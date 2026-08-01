#!/bin/bash
# Build the mzPeak C++ reader (writer_test) on spock. It has a Rust core, so this
# installs a Rust toolchain if absent. Everything under /scratch (NEVER /home:
# quota). Gating dependency for the streaming extractor.
set -uo pipefail
export CARGO_HOME=/scratch/kohlbach/home/.cargo
export RUSTUP_HOME=/scratch/kohlbach/home/.rustup
export PATH=$CARGO_HOME/bin:/scratch/kohlbach/mamba/envs/odia/bin:$PATH
M=/scratch/kohlbach/opendialyzer/ext/mzpeak
cd "$M"

echo "[$(date +%T)] rust: $(command -v cargo || echo absent)"
if ! command -v cargo >/dev/null; then
  echo "[$(date +%T)] installing rustup (offline-safe dirs on /scratch) ..."
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --no-modify-path --profile minimal 2>&1 | tail -5
fi
echo "[$(date +%T)] cargo=$(cargo --version 2>&1) meson=$(meson --version 2>&1)"

echo "[$(date +%T)] meson setup ..."
rm -rf build
meson setup build 2>&1 | tail -8
echo "[$(date +%T)] meson compile ..."
meson compile -C build 2>&1 | tail -25

echo "[$(date +%T)] artifacts:"
find build -name "*.a" -o -name "*.so" 2>/dev/null | head
echo "[$(date +%T)] DONE"
