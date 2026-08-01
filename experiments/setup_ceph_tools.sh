#!/usr/bin/env bash
# Install shared tools to /ceph/ibmi/abi/oliver so every IBMI node can use them
# (/scratch is node-local -- that is how one node became a hard dependency before).
set -uo pipefail
export HOME=/scratch/$USER/home TMPDIR=/scratch/$USER/tmp
B=/ceph/ibmi/abi/oliver
mkdir -p "$B/opt" "$B/bin" "$B/src" "$TMPDIR"
log(){ printf '[%s] %s\n' "$(date +%T)" "$*"; }

# ---- DIA-NN (open-source build from GitHub releases) -------------------------
log "resolving latest DIA-NN release"
API=https://api.github.com/repos/vdemichev/DiaNN/releases
URL=$(curl -sS "$API" 2>/dev/null | grep -oE '"browser_download_url": *"[^"]*"' | cut -d'"' -f4 \
      | grep -iE 'linux|\.zip$|\.tar\.gz$' | head -1)
if [[ -z ${URL:-} ]]; then
  log "no asset via API; listing what IS offered:"
  curl -sS "$API" 2>/dev/null | grep -oE '"browser_download_url": *"[^"]*"' | cut -d'"' -f4 | head -8
else
  log "DIA-NN asset: $URL"
  F=$B/src/$(basename "$URL")
  [[ -s $F ]] || curl -sSL --retry 5 -o "$F" "$URL"
  mkdir -p "$B/opt/diann"
  case "$F" in
    *.zip)     unzip -qo "$F" -d "$B/opt/diann" ;;
    *.tar.gz)  tar -xzf "$F" -C "$B/opt/diann" ;;
  esac
  BIN=$(find "$B/opt/diann" -maxdepth 3 -type f \( -name 'diann*' -o -name 'DiaNN*' \) -perm -u+x 2>/dev/null | head -1)
  [[ -n ${BIN:-} ]] && { chmod +x "$BIN"; ln -sf "$BIN" "$B/bin/diann"; log "diann -> $BIN"; }
fi

# ---- mzPeakConverter ---------------------------------------------------------
log "installing mzpeak-convert"
if [[ -x /scratch/$USER/src/mzPeakConverter/target/release/mzpeak-convert ]]; then
  cp /scratch/$USER/src/mzPeakConverter/target/release/mzpeak-convert "$B/bin/"
  log "mzpeak-convert installed from the existing build"
else
  log "no local build to copy"
fi

log "contents:"; ls -la "$B/bin" 2>/dev/null | tail -5
