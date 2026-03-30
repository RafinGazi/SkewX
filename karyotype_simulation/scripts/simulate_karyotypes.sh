#!/bin/bash
set -euo pipefail
module load SAMtools

BASE_DIR="/mnt/Genomics/Lab/HEAL/X_chr/Rafin/SkewX/karyotype_simulation"
INT_DIR="${BASE_DIR}/intermediate"
OUT_DIR="${BASE_DIR}/output"
mkdir -p "$OUT_DIR"

echo "Starting karyotype simulation..."

AUTO="${INT_DIR}/GM19312_auto.bam"
X="${INT_DIR}/GM19312_chrX.bam"
Y="${INT_DIR}/GM19312_chrY.bam"

############################################
# Helper: rename read names with a suffix
# Usage: rename_reads INPUT OUTPUT SUFFIX
############################################
rename_reads () {
    INPUT=$1
    OUTPUT=$2
    SUFFIX=$3
    echo "  Renaming reads in $(basename $INPUT) with suffix _${SUFFIX}..."
    samtools view -h "$INPUT" | \
        awk -v suf="_${SUFFIX}" 'BEGIN{OFS="\t"} /^@/{print; next} {$1=$1 suf; print}' | \
        samtools view -b -o "$OUTPUT"
    samtools index "$OUTPUT"
}

############################################
# Helper: merge, sort, index, cleanup
############################################
build_and_index () {
    NAME=$1
    shift
    echo "Building $NAME..."
    samtools merge -f "${OUT_DIR}/${NAME}.bam" "$@"
    samtools sort "${OUT_DIR}/${NAME}.bam" -o "${OUT_DIR}/${NAME}_sorted.bam"
    samtools index "${OUT_DIR}/${NAME}_sorted.bam"
    rm "${OUT_DIR}/${NAME}.bam"
    echo "$NAME done."
}

############################################
# Pre-build renamed copies to avoid
# read name collisions when merging same BAM
############################################
echo "Creating renamed copies of chrX and chrY..."
rename_reads "$X" "${OUT_DIR}/chrX_copy2.bam" "copy2"
rename_reads "$X" "${OUT_DIR}/chrX_copy3.bam" "copy3"
rename_reads "$Y" "${OUT_DIR}/chrY_copy2.bam" "copy2"

############################################
# Simulations
############################################

# XY (baseline — no duplicates, no renaming needed)
build_and_index "XY" "$AUTO" "$X" "$Y"

# XX — 2x chrX (use original + renamed copy)
build_and_index "XX" "$AUTO" "$X" "${OUT_DIR}/chrX_copy2.bam"

# XO — half coverage chrX, no Y
samtools view -b -s 0.5 "$X" -o "${OUT_DIR}/X_half.bam"
samtools index "${OUT_DIR}/X_half.bam"
build_and_index "XO" "$AUTO" "${OUT_DIR}/X_half.bam"

# XXX — 3x chrX
build_and_index "XXX" "$AUTO" "$X" "${OUT_DIR}/chrX_copy2.bam" "${OUT_DIR}/chrX_copy3.bam"

# XXY — 2x chrX + 1x chrY
build_and_index "XXY" "$AUTO" "$X" "${OUT_DIR}/chrX_copy2.bam" "$Y"

# XYY — 1x chrX + 2x chrY
build_and_index "XYY" "$AUTO" "$X" "$Y" "${OUT_DIR}/chrY_copy2.bam"

############################################
# Cleanup renamed intermediates
############################################
echo "Cleaning up intermediate renamed BAMs..."
rm -f "${OUT_DIR}/chrX_copy2.bam" "${OUT_DIR}/chrX_copy2.bam.bai"
rm -f "${OUT_DIR}/chrX_copy3.bam" "${OUT_DIR}/chrX_copy3.bam.bai"
rm -f "${OUT_DIR}/chrY_copy2.bam" "${OUT_DIR}/chrY_copy2.bam.bai"
rm -f "${OUT_DIR}/X_half.bam"     "${OUT_DIR}/X_half.bam.bai"

echo "Simulation complete."