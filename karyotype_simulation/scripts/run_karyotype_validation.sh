#!/bin/bash
set -euo pipefail

# =========================

# LOAD MODULES

# =========================

module load samtools || true

# =========================

# PATH SETUP

# =========================

BASE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$BASE/output"
MOS="$BASE/mosdepth"
RES="$BASE/results"

mkdir -p "$MOS" "$RES"

echo "==== STEP 0: INPUT CHECK ===="
echo "BASE=$BASE"
echo "OUT=$OUT"

if ! ls "$OUT"/*_sorted.bam 1> /dev/null 2>&1; then
echo "ERROR: No BAM files found in $OUT"
exit 1
fi

ls -lh "$OUT"

shopt -s nullglob

# =========================

# STEP 1: MOSDEPTH

# =========================

echo "==== STEP 1: MOSDEPTH ===="

for BAM in "$OUT"/*_sorted.bam; do
NAME=$(basename "$BAM" .bam)
echo "[mosdepth] Processing $NAME"

```
singularity exec -B "$BASE" \
    docker://quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0 \
    mosdepth -t 4 -b 1000000 "$MOS/$NAME" "$BAM"
```

done

# =========================

# STEP 2: KARYOTYPE INFERENCE

# =========================

echo "==== STEP 2: KARYOTYPE INFERENCE ===="

if ! ls "$MOS"/*.mosdepth.summary.txt 1> /dev/null 2>&1; then
echo "ERROR: No mosdepth outputs found"
exit 1
fi

cd "$RES"

for SUMMARY in "$MOS"/*.mosdepth.summary.txt; do
NAME=$(basename "$SUMMARY" .mosdepth.summary.txt)
PREFIX="$MOS/$NAME"

```
echo "[infer] Processing $NAME"

if python3 "/mnt/Genomics/Lab/HEAL/X_chr/Rafin/SkewX/bin/infer_karyotype_v2.py" \
    "$PREFIX.regions.bed.gz" \
    "$SUMMARY" \
    "$NAME" \
    "${NAME}_karyotype.tsv"; then

    echo "[INFO] $NAME passed (XX)"

else
    EXIT_CODE=$?

    if [ $EXIT_CODE -eq 2 ]; then
        echo "[INFO] $NAME skipped (non-XX karyotype — expected in validation)"
    else
        echo "[ERROR] $NAME failed during karyotype inference"
        exit 1
    fi
fi
```

done

# =========================

# STEP 3: GENERATE HTML REPORT

# =========================

echo "==== STEP 3: GENERATE HTML REPORT ===="

cd "$BASE/scripts"
python3 report.py

echo "==== ALL DONE ===="
echo "Results directory: $RES"
echo "Open report: $RES/report.html"
