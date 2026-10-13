#!/bin/bash
#SBATCH --job-name=magus_final
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --array=1-2124%10

# Align final locus sequences with MAGUS.
# Submit this script from the project working directory.

LOCI_NAME_FILE="2-loci_filter/loci_name_filter.log"
INPUT_DIR="1-raw_loci_final/faa"
OUTPUT_DIR="3-faa_aln_final"
FAILED_DIR="logs/magus_failed"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTPUT_DIR" "$FAILED_DIR"

LOCI_FILTER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LOCI_NAME_FILE")
INPUT_FILE="${INPUT_DIR}/${LOCI_FILTER}.faa"
OUTPUT_FILE="${OUTPUT_DIR}/${LOCI_FILTER}.faa"
WORKDIR="${OUTPUT_DIR}/magus_working_dir_${LOCI_FILTER}"

echo "Starting $LOCI_FILTER (task $SLURM_ARRAY_TASK_ID) at $(date)"

rm -rf "$WORKDIR"

if [[ -s "$INPUT_FILE" ]]; then
    magus \
        -i "$INPUT_FILE" \
        -o "$OUTPUT_FILE" \
        -np "$THREADS" \
        -d "$WORKDIR"

    if [[ -s "$OUTPUT_FILE" ]]; then
        echo "Locus $LOCI_FILTER aligned successfully"
        rm -rf "$WORKDIR"
    else
        echo "Locus $LOCI_FILTER alignment failed" \
            | tee -a "${FAILED_DIR}/list.txt"

        if [[ -d "$WORKDIR" ]]; then
            mv "$WORKDIR" \
               "${FAILED_DIR}/${LOCI_FILTER}_working_dir_$(date +%s)"
        fi
    fi
else
    echo "Input file missing: $INPUT_FILE" \
        | tee -a "${FAILED_DIR}/list.txt"
fi
