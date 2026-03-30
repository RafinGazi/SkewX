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
# Helper: merge, sort, index, cleanup
############################################
build_and_index () {
    NAME=$1
    shift
    echo "Building $NAME..."
    samtools merge "${OUT_DIR}/${NAME}.bam" "$@"
    samtools sort "${OUT_DIR}/${NAME}.bam" -o "${OUT_DIR}/${NAME}_sorted.bam"
    samtools index "${OUT_DIR}/${NAME}_sorted.bam"
    rm "${OUT_DIR}/${NAME}.bam"
    echo "$NAME done."
}

############################################
# Simulations
############################################

# XY — 1x chrX + 1x chrY
build_and_index "XY" "$AUTO" "$X" "$Y"

# XX — 2x chrX
build_and_index "XX" "$AUTO" "$X" "$X"

# XO — half coverage chrX, no Y
samtools view -b -s 0.5 "$X" -o "${OUT_DIR}/X_half.bam"
samtools index "${OUT_DIR}/X_half.bam"
build_and_index "XO" "$AUTO" "${OUT_DIR}/X_half.bam"
rm -f "${OUT_DIR}/X_half.bam" "${OUT_DIR}/X_half.bam.bai"

# XXX — 3x chrX
build_and_index "XXX" "$AUTO" "$X" "$X" "$X"

# XXY — 2x chrX + 1x chrY
build_and_index "XXY" "$AUTO" "$X" "$X" "$Y"

# XYY — 1x chrX + 2x chrY
build_and_index "XYY" "$AUTO" "$X" "$Y" "$Y"

echo "Simulation complete."