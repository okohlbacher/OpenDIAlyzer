#!/bin/bash
# Rebuild the agxt library with the C2/C5 fixes: timsTOF instrument index (2),
# NCE 35 (measured from analysis.tdf: 25.8-46.8 eV, mean 35.4), and protein-N-term
# Acetyl gated to actual protein N-terminal peptides.
# raw -> m/z 300-1200 trim -> AssayGenerator -> DecoyGenerator.
set -uo pipefail
O=/scratch/kohlbach/opendialyzer/bench/agxt_fixed
MD=/scratch/kohlbach/opendialyzer/build/openms-onnx/src/openms/share/OpenMS/models
OMS=/home/kohlbach/openms3/bin
mkdir -p $O

echo "[$(date +%T)] raw generation (timsTOF/NCE35, acetyl gated) ..."
( export LD_LIBRARY_PATH=/scratch/kohlbach/opendialyzer/build/openms-onnx/lib
  cd /scratch/kohlbach/opendialyzer/OpenDIAlyzer/build-s1
  /usr/bin/time -v ./OpenDIALibGen -fasta /scratch/agxt/diann/agxt_variants.fasta \
    -out $O/raw.tsv -model_dir $MD -threads 128 -max_charge 4 -max_var_mods 1 -raw ) \
  > $O/raw_gen.log 2>&1 || { echo "FATAL raw gen"; exit 1; }
echo "[$(date +%T)] raw: $(du -h $O/raw.tsv | cut -f1)"

# Match DIA-NN's search space: precursor m/z 300-1200 (its lib-gen setting).
echo "[$(date +%T)] trim to m/z 300-1200 ..."
awk -F'\t' -v OFS='\t' 'NR==1{print;next}{mz=$1+0; if(mz>=300&&mz<=1200) print}' \
  $O/raw.tsv > $O/raw_trim.tsv
echo "[$(date +%T)] trimmed: $(du -h $O/raw_trim.tsv | cut -f1)"

export LD_LIBRARY_PATH=/home/kohlbach/openms3/lib
echo "[$(date +%T)] assay generator (P1: 6-12 transitions, product floor 200) ..."
# P1 fragment-selection fix applied at the TOPP level (the path used at scale;
# the in-process refine_and_decoy is bypassed here). Defaults 350/6 gut a DIA
# library -- floor 150 keeps low-mass base peaks, 6-12 matches DIA-NN's count.
$OMS/OpenSwathAssayGenerator -in $O/raw_trim.tsv -out $O/assay.tsv \
  -min_transitions 6 -max_transitions 12 \
  -product_lower_mz_limit 200 -product_upper_mz_limit 2000 \
  -precursor_lower_mz_limit 300 -precursor_upper_mz_limit 1200 > $O/assay.log 2>&1 \
  || { echo "FATAL assay"; exit 2; }
echo "[$(date +%T)] decoy generator ..."
$OMS/OpenSwathDecoyGenerator -in $O/assay.tsv -out $O/library.tsv \
  -method shuffle -switchKR true -min_decoy_fraction 0.1 > $O/decoy.log 2>&1 \
  || { echo "FATAL decoy"; exit 3; }
echo "[$(date +%T)] LIBRARY READY: $(du -h $O/library.tsv | cut -f1)"
awk -F'\t' 'NR>1{print $25}' $O/library.tsv | sort | uniq -c
