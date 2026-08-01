#!/bin/bash
# Regenerate the OpenSWATH patch from ext/OpenMS, and FAIL if it changed.
#
# ext/ is gitignored, so an edit to the OpenMS core leaves no trace in this repository unless the
# patch is regenerated. Running this in CI (or before a commit) turns "someone forgot" from a silent
# loss into a failed check. PEPTDEEP is excluded: it has its own tracked patch.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OUT="$ROOT/vendored-patches/OpenMS/opendialyzer-openswath.patch"
[ -d "$ROOT/ext/OpenMS" ] || { echo "ext/OpenMS absent — nothing to check"; exit 0; }

TMP=$(mktemp); trap 'rm -f "$TMP"' EXIT
git -C "$ROOT/ext/OpenMS" diff -- . ':(exclude)*PEPTDEEP*' > "$TMP"

if [ "${1:-}" = "--check" ]; then
  if ! diff -q "$TMP" "$OUT" >/dev/null 2>&1; then
    echo "DRIFT: ext/OpenMS does not match $(basename "$OUT")." >&2
    diff -u "$OUT" "$TMP" | head -40 >&2
    echo >&2; echo "Run $0 to update it, and document the change in vendored-patches/README.md." >&2
    exit 1
  fi
  echo "patch is in sync with ext/OpenMS"
else
  cp "$TMP" "$OUT"
  echo "regenerated $OUT"
  git -C "$ROOT" diff --stat -- "$OUT" || true
fi
