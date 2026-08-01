#!/usr/bin/env bash
# PXD017703 (Meier diaPASEF). Smallest useful arm FIRST so the three-way can start sooner:
# HeLa_Evosep is 9 short-gradient .d runs at ~1.2 GB each (11.08 GB total) versus the
# TwoProteome zip's 68.92 GB. Both ship a prebuilt PQP library and the published pyprophet
# result, so the ~650 GB of ddaPASEF library raws are NOT needed.
set -uo pipefail
export HOME=/scratch/$USER/home TMPDIR=/scratch/$USER/tmp
D=/ceph/ibmi/abi/data/dia_reference/PXD017703_diaPASEF
mkdir -p "$D"
R=ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/12/PXD017703   # path from the v3 API
log(){ printf '[%s] %s\n' "$(date +%T)" "$*"; }
# smallest first, then the big ground-truth zip
for f in HeLa_Evosep_pqp_library.zip HeLa_Evosep_pyprophet_export.zip HeLa_Evosep_diaPASEF_RAW.zip \
         TwoProteome_pqp_library.zip TwoProteome_pyprophet_export.zip TwoProteome_diaPASEF_raw.zip; do
  t=$D/$f
  # -C - resumes; PRIDE drops the connection repeatedly on the multi-GB zips.
  if curl -sS -L --retry 20 --retry-delay 20 --retry-all-errors -C - -o "$t" "$R/$f"; then
    log "ok $f ($(du -h "$t" | cut -f1))"
  else
    log "FAILED $f (partial: $(du -h "$t" 2>/dev/null | cut -f1))"
  fi
done
log "total $(du -sh "$D" | cut -f1)"
