#!/usr/bin/env bash
# PRIDE over Aspera (fasp). HTTPS/FTP keep dropping the multi-GB zips with curl (18);
# fasp is EBI's supported path for large transfers and resumes properly.
set -uo pipefail
export HOME=/scratch/$USER/home TMPDIR=/scratch/$USER/tmp
B=/ceph/ibmi/abi/oliver
D=/ceph/ibmi/abi/data/dia_reference/PXD017703_diaPASEF
mkdir -p "$D" "$TMPDIR"
log(){ printf '[%s] %s\n' "$(date +%T)" "$*"; }

ASCP=$(command -v ascp || true)
if [[ -z $ASCP ]]; then
  log "installing aspera-cli into $B/envs/aspera"
  /scratch/$USER/micromamba/micromamba create -y -p "$B/envs/aspera" -c conda-forge -c hcc aspera-cli \
    > "$TMPDIR/aspera_install.log" 2>&1
  ASCP=$(find "$B/envs/aspera" -name ascp -type f -perm -u+x 2>/dev/null | head -1)
fi
[[ -z ${ASCP:-} ]] && { log "FATAL: no ascp available"; tail -5 "$TMPDIR/aspera_install.log"; exit 1; }
log "ascp: $ASCP"

KEY=$(find "$(dirname "$(dirname "$ASCP")")" -name "asperaweb_id_dsa.openssh" 2>/dev/null | head -1)
[[ -z ${KEY:-} ]] && KEY=$(find "$B/envs/aspera" -name "asperaweb_id_dsa.openssh" 2>/dev/null | head -1)
log "key: ${KEY:-<none found>}"

# -QT unlimited-rate fair transfer, -l 500m cap, -k2 resume by sparse checksum, -P33001 fasp port
for f in HeLa_Evosep_diaPASEF_RAW.zip TwoProteome_diaPASEF_raw.zip; do
  log "GET $f"
  "$ASCP" -QT -l 500m -P33001 -k2 ${KEY:+-i "$KEY"} \
    "prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2020/12/PXD017703/$f" "$D/" \
    && log "ok $f ($(du -h "$D/$f" | cut -f1))" || log "FAILED $f"
done
log "total $(du -sh "$D" | cut -f1)"
