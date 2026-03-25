#!/bin/bash

set -euo pipefail

module load SAMtools

BASE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$BASE/output"
MOS="$BASE/mosdepth"
RES="$BASE/results"

mkdir -p "$MOS" "$RES"

echo "==== STEP 1: MOSDEPTH ===="

echo "BASE=$BASE"
echo "OUT=$OUT"
ls -lh "$OUT"

for BAM in $OUT/*_sorted.bam; do
    NAME=$(basename "$BAM" .bam)
    echo "Processing $NAME"

    singularity exec -B /mnt/Genomics/Lab/HEAL/X_chr/Rafin \
    	docker://quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0 \
    	mosdepth -t 4 -b 1000000 "$MOS/$NAME" "$BAM"
done

echo "==== STEP 2: KARYOTYPE INFERENCE ===="

for SUMMARY in $MOS/*.mosdepth.summary.txt; do
    NAME=$(basename "$SUMMARY" .mosdepth.summary.txt)
    PREFIX="$MOS/$NAME"

    python3 /mnt/Genomics/Lab/HEAL/X_chr/Rafin/SkewX/bin/infer_karyotype_v2.py \
        "$PREFIX.regions.bed.gz" \
        "$SUMMARY" \
        "$NAME" \
        "$RES/${NAME}_karyotype.tsv"
done

echo "==== STEP 3: COHORT QC ===="

python3 /mnt/Genomics/Lab/HEAL/X_chr/Rafin/SkewX/bin/cohort_karyotype_qc.py \
    $RES/*.tsv \
    --out_prefix "$RES/cohort"

echo "==== STEP 4: GENERATE HTML REPORTS ===="

cd "$BASE/scripts"
python3 report.py

echo "==== ALL DONE ===="
