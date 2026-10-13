# BUSCO protein alignment and alignment trimming
# Run from the BUSCO extraction working directory.
# Required software: MAGUS, trimAl

THREADS=8
export OMP_NUM_THREADS=$THREADS

# Align protein sequences with MAGUS
mkdir -p 3-faa_align logs/magus_failed
mapfile -t LOCI_NAME_FILTER < 2-loci_filter/loci_name_filter.log
: > logs/magus_failed/list.txt

for LOCI_FILTER in "${LOCI_NAME_FILTER[@]}"
do
    WORKDIR="3-faa_align/magus_working_dir_${LOCI_FILTER}"
    rm -rf "$WORKDIR"

    magus \
        -i "2-loci_filter/faa/${LOCI_FILTER}.faa" \
        -o "3-faa_align/${LOCI_FILTER}.faa" \
        -np "$THREADS" \
        -d "$WORKDIR"

    if test -s "3-faa_align/${LOCI_FILTER}.faa"; then
        echo "Locus $LOCI_FILTER aligned successfully"
        rm -rf "$WORKDIR"
    else
        echo "$LOCI_FILTER" >> logs/magus_failed/list.txt
        if [ -d "$WORKDIR" ]; then
            mv "$WORKDIR" "logs/magus_failed/${LOCI_FILTER}_working_dir"
        fi
    fi
 done

# Trim protein alignments and back-translate to nucleotide alignments
mkdir -p 4-loci_trim/fna 4-loci_trim/faa

for LOCI_FILTER in "${LOCI_NAME_FILTER[@]}"
do
    [ -s "3-faa_align/${LOCI_FILTER}.faa" ] || continue

    trimal \
        -in "3-faa_align/${LOCI_FILTER}.faa" \
        -out "4-loci_trim/faa/${LOCI_FILTER}.aa.fas" \
        -automated1

    trimal \
        -in "3-faa_align/${LOCI_FILTER}.faa" \
        -out "4-loci_trim/fna/${LOCI_FILTER}.nuc.fas" \
        -automated1 \
        -backtrans "2-loci_filter/fna/${LOCI_FILTER}.fna"
 done
