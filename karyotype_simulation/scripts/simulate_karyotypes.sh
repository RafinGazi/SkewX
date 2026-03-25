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

build_and_index () {
    NAME=$1
    shift

    echo "Building $NAME..."

    samtools merge "${OUT_DIR}/${NAME}.bam" "$@"
    samtools sort "${OUT_DIR}/${NAME}.bam" -o "${OUT_DIR}/${NAME}_sorted.bam"
    samtools index "${OUT_DIR}/${NAME}_sorted.bam"

    rm "${OUT_DIR}/${NAME}.bam"
}

############################################
# Simulations
############################################

# XY (baseline)
build_and_index "XY" "$AUTO" "$X" "$Y"

# XX
build_and_index "XX" "$AUTO" "$X" "$X"

# XO
samtools view -b -s 0.5 "$X" -o "${OUT_DIR}/X_half.bam"
samtools index "${OUT_DIR}/X_half.bam"
build_and_index "XO" "$AUTO" "${OUT_DIR}/X_half.bam"

# XXX
build_and_index "XXX" "$AUTO" "$X" "$X" "$X"

# XXY
build_and_index "XXY" "$AUTO" "$X" "$X" "$Y"

# XYY
build_and_index "XYY" "$AUTO" "$X" "$Y" "$Y"

echo "Simulation complete."