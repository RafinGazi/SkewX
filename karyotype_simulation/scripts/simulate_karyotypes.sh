#!/bin/bash

set -euo pipefail

module load SAMtools

BASE_DIR="/mnt/Genomics/Lab/HEAL/X_chr/Rafin/karyotype_simulation"
INT_DIR="${BASE_DIR}/intermediate"
OUT_DIR="${BASE_DIR}/output"

mkdir -p "$OUT_DIR"

echo "Starting karyotype simulation..."

############################################
# Base components
############################################

AUTO="${INT_DIR}/GM19462_auto.bam"
X1="${INT_DIR}/GM19462_chrX.bam"
X2="${INT_DIR}/HG01363_chrX.bam"
Y1="${INT_DIR}/GM19312_chrY.bam"
Y2="${INT_DIR}/GM18865_chrY.bam"

############################################
# Helper function
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

# XX (baseline)
build_and_index "XX" "$AUTO" "$X1" "$X2"

# XY (baseline)
build_and_index "XY" "$AUTO" "$X1" "$Y1"

# XO (half X)
samtools view -b -s 0.5 "$X1" -o "${OUT_DIR}/X_half.bam"
samtools index "${OUT_DIR}/X_half.bam"

build_and_index "XO" "$AUTO" "${OUT_DIR}/X_half.bam"

# XXX
build_and_index "XXX" "$AUTO" "$X1" "$X2" "$X2"

# XXY
build_and_index "XXY" "$AUTO" "$X1" "$X2" "$Y1"

# XYY
build_and_index "XYY" "$AUTO" "$X1" "$Y1" "$Y2"

echo "All karyotypes simulated successfully."
