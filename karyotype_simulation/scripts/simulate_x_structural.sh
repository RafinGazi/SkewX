#!/bin/bash

set -euo pipefail

module load SAMtools

BASE_DIR="/mnt/Genomics/Lab/HEAL/X_chr/Rafin/SkewX/karyotype_simulation"
INT_DIR="${BASE_DIR}/intermediate"
OUT_DIR="${BASE_DIR}/output"

mkdir -p "$OUT_DIR"

echo "Starting X structural variant simulations..."

############################################
# Base files
############################################

AUTO="${INT_DIR}/GM19462_auto.bam"
X="${INT_DIR}/GM19462_chrX.bam"

############################################
# Define regions
############################################

XP_REGION="chrX:1-60000000"
XQ_REGION="chrX:60000001-156040895"

############################################
# Helper function
############################################

build_and_index () {
    NAME=$1
    shift

    echo "Building $NAME..."

    # Merge WITHOUT sorting
    samtools merge "${OUT_DIR}/${NAME}_sorted.bam" "$@"

    # Try indexing
    if samtools index "${OUT_DIR}/${NAME}_sorted.bam"; then
        echo "Index successful for $NAME (sorting not required)"
    else
        echo "Index failed → sorting required for $NAME"

        # Fallback: sort then index
        samtools sort "${OUT_DIR}/${NAME}_sorted.bam" -o "${OUT_DIR}/${NAME}_tmp.bam"
        mv "${OUT_DIR}/${NAME}_tmp.bam" "${OUT_DIR}/${NAME}_sorted.bam"
        samtools index "${OUT_DIR}/${NAME}_sorted.bam"
    fi

    # Quick integrity check
    samtools quickcheck "${OUT_DIR}/${NAME}_sorted.bam" || {
        echo "ERROR: BAM failed integrity check for $NAME"
        exit 1
    }

    echo "$NAME done."
}

############################################
# 1. Xp deletion (remove short arm)
############################################

samtools view -b "$X" $XQ_REGION -o "${OUT_DIR}/Xq_only.bam"
samtools index "${OUT_DIR}/Xq_only.bam"

build_and_index "Xp_deleted" "$AUTO" "${OUT_DIR}/Xq_only.bam"

############################################
# 2. Xq deletion (remove long arm)
############################################

samtools view -b "$X" $XP_REGION -o "${OUT_DIR}/Xp_only.bam"
samtools index "${OUT_DIR}/Xp_only.bam"

build_and_index "Xq_deleted" "$AUTO" "${OUT_DIR}/Xp_only.bam"

############################################
# 3. Partial deletion (50% drop in Xp)
############################################

# Xp half coverage
samtools view -b "$X" $XP_REGION | samtools view -b -s 0.5 -o "${OUT_DIR}/Xp_half.bam"
samtools index "${OUT_DIR}/Xp_half.bam"

# Xq full
samtools view -b "$X" $XQ_REGION -o "${OUT_DIR}/Xq_full.bam"
samtools index "${OUT_DIR}/Xq_full.bam"

build_and_index "X_partial" "$AUTO" "${OUT_DIR}/Xp_half.bam" "${OUT_DIR}/Xq_full.bam"

############################################
# Cleanup intermediates (optional)
############################################

rm -f "${OUT_DIR}/Xq_only.bam" "${OUT_DIR}/Xq_only.bam.bai"
rm -f "${OUT_DIR}/Xp_only.bam" "${OUT_DIR}/Xp_only.bam.bai"
rm -f "${OUT_DIR}/Xp_half.bam" "${OUT_DIR}/Xp_half.bam.bai"
rm -f "${OUT_DIR}/Xq_full.bam" "${OUT_DIR}/Xq_full.bam.bai"

echo "All structural simulations completed."
