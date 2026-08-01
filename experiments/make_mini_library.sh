#!/bin/bash
# Build a small benchmark library by keeping every Nth transition group, PLUS
# every group whose peptide is a CiRT anchor. The anchors are not optional: a
# naive 1/N sample drops them and OpenSWATH then dies with "iRT calibration
# failed: insufficient RT coverage after outlier removal".
# Usage: make_mini_library.sh <full_library.tsv> <out.tsv> [N]
set -uo pipefail
L=$1
OUT=$2
N=${3:-240}
CIRT=/home/kohlbach/openms3/share/OpenMS/CHEMISTRY/cirtkit.tsv
TMP=$(dirname "$OUT")/cirt.seq
awk -F'\t' 'NR>1{print $7}' "$CIRT" | sort -u > "$TMP"

awk -F'\t' -v OFS='\t' -v N="$N" '
  NR==FNR { cirt[$1]=1; next }
  FNR==1  { print; next }
  {
    g=$23
    if (!(g in seen)) seen[g] = ((++n % N)==0) || ($7 in cirt)
    if (seen[g]) print
  }
' "$TMP" "$L" > "$OUT"

echo "rows:   $(wc -l < "$OUT")"
echo "groups: $(awk -F'\t' 'NR>1{print $23}' "$OUT" | sort -u | wc -l)"
echo "targets/decoys:"; awk -F'\t' 'NR>1{print $25}' "$OUT" | sort | uniq -c
echo "CiRT peptides: $(awk -F'\t' 'NR>1{print $7}' "$OUT" | sort -u | grep -Fxf "$TMP" | wc -l) / $(wc -l < "$TMP")"
awk -F'\t' 'NR>1{r=$6+0; if(r<mn||NR==2)mn=r; if(r>mx)mx=r} END{printf "RT span: %.4f .. %.4f\n", mn, mx}' "$OUT"
