# 4.2 Supermatrix concatenation
# Example shown for SRH_taxa75. The same procedure was applied to other datasets.
# Internal labels taxa50/75/90 correspond to matrix50/75/90 in the final manuscript.

DATASET="SRH_taxa75"
ALIGN_DIR="../busco_extraction/4-loci_trim_final.clean/faa"
LOCI_LIST="loci_list/loci_list_${DATASET}.txt"
FASCONCAT="FASconCAT-G_v1.05.1.pl"

mkdir -p "$DATASET"

for LOCUS in $(cat "$LOCI_LIST"); do
    cp "${ALIGN_DIR}/${LOCUS}.aa.fas" "${DATASET}/"
done

(
    cd "$DATASET" || exit 1
    perl "../${FASCONCAT}" -s -l -p
)

# Create the secondary dataset retaining sequences with >=70% completeness.

./4.2-script_filter_by_completeness.sh \
    -i "$DATASET" \
    -o "${DATASET}/fas_len70" \
    -l "$LOCI_LIST" \
    -c 70 \
    -r "$DATASET"

(
    cd "${DATASET}/fas_len70" || exit 1
    perl "../../${FASCONCAT}" -s -l -p
)
