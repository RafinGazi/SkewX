#!/bin/bash

set -euo pipefail

DATA_DIR="/mnt/Genomics/Lab/HEAL/X_chr/Rafin/data"
OUT_DIR="/mnt/Genomics/Lab/HEAL/X_chr/Rafin/karyotype_simulation/intermediate"

mkdir -p "$OUT_DIR"

SAMPLE="GM19312"   # male sample

echo "Processing $SAMPLE..."

INPUT_BAM="${DATA_DIR}/${SAMPLE}_chrX_Y_18_21.bam"

# Output files
X_BAM="${OUT_DIR}/${SAMPLE}_chrX.bam"
Y_BAM="${OUT_DIR}/${SAMPLE}_chrY.bam"
AUTO_BAM="${OUT_DIR}/${SAMPLE}_auto.bam"

# Extract chrX
samtools view -b "$INPUT_BAM" chrX -o "$X_BAM"

# Extract chrY
samtools view -b "$INPUT_BAM" chrY -o "$Y_BAM"

# Extract autosomes (chr18 + chr21)
samtools view -b "$INPUT_BAM" chr18 chr21 -o "$AUTO_BAM"

# Index all
samtools index "$X_BAM"
samtools index "$Y_BAM"
samtools index "$AUTO_BAM"

echo "Extraction complete."