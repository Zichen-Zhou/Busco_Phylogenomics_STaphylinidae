# 3.2 Non-ortholog detection and removal
# Run these blocks from the tree-based filtering working directory.

LOCI_LIST_FILE="../2-loci_filter/loci_name_filter.log"
APPLY_SCRIPT="3.2-script_apply_treeshrink_removal_4.0.py"


# 3.2.1 Preliminary long-branch detection with TreeShrink

echo "Locus,SeqID" > treeshrink_rmseq.csv

MODE="per-gene"
QUANTILES="0.4"
TREE_DIR="genetrees_ori"
LOCI_LIST=$(cat "$LOCI_LIST_FILE")

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"
    TMP_DIR="tmp_treeshrink/${LOCI}"

    if [[ -s "$TREE" ]]; then
        mkdir -p "$TMP_DIR"
        cp "$TREE" "$TMP_DIR/input.tree"

        (
            cd "$TMP_DIR" || exit 1
            run_treeshrink.py -t input.tree -m "$MODE" -q "$QUANTILES" -O output > treeshrink.log 2>&1
        )

        RMSEQ_FILE="${TMP_DIR}/input_treeshrink/output.txt"
        if [[ -f "$RMSEQ_FILE" ]]; then
            while IFS=$'\t' read -r -a SEQS; do
                for ID in "${SEQS[@]}"; do
                    if [[ -n "$ID" ]]; then
                        echo "${LOCI},${ID}" >> treeshrink_rmseq.csv
                    fi
                done
            done < "$RMSEQ_FILE"
        else
            echo "[Warning] TreeShrink output not found for $LOCI" >&2
        fi
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done


# 3.2.2 Apply preliminary TreeShrink removal

python3 "$APPLY_SCRIPT" \
    -l "$LOCI_LIST_FILE" \
    -r treeshrink_rmseq.csv \
    -i ../1-raw_loci \
    -o ../1-raw_loci_rmTS \
    -c sequence_number_rmTS.log \
    -n loci_name_rmTS.log

# Gene trees were rebuilt after this removal using the procedure in
# 3.1-code_genetree_building.sh.


# 3.2.3 Average branch support (ABS) for preliminary checking

LOCI_LIST=$(cat "$LOCI_LIST_FILE")
TREE_DIR="genetrees_rmTS"
OUT_CSV="FilterPara_rmTS/tree_based_stat_rmTS.csv"

mkdir -p FilterPara_rmTS
echo "Locus,ABS" > "$OUT_CSV"

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"

    if [[ -s "$TREE" ]]; then
        ABS=$(sed 's/)/\n/g' "$TREE" \
            | cut -d ":" -f 1 \
            | sed '1d' \
            | sed 's/;//g' \
            | awk NF \
            | awk '{a+=$1} END {if (NR>0) print a/NR; else print "NA"}')
        echo "${LOCI},${ABS}" >> "$OUT_CSV"
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done

csvtk plot hist -f ABS --format pdf \
    -o FilterPara_rmTS/03_ABS_rmTS.pdf \
    < FilterPara_rmTS/tree_based_stat_rmTS.csv


# 3.2.4 Treeness for non-ortholog checking

TREE_DIR="genetrees_rmTS"
TMP_DIR="FilterPara_rmTS/tmp5.0.3_treeness"
OUT_CSV="${TMP_DIR}/treeness.csv"

mkdir -p "$TMP_DIR"
echo "Locus,treeness" > "$OUT_CSV"

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"

    if [[ -s "$TREE" ]]; then
        treeness=$(phykit treeness "$TREE" 2>/dev/null)
        echo "${LOCI},${treeness}" >> "$OUT_CSV"
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done

csvtk join -H -f Locus \
    FilterPara_rmTS/tree_based_stat_rmTS.csv \
    "$OUT_CSV" \
    > FilterPara_rmTS/tree_based_stat_rmTS.updated.csv

mv FilterPara_rmTS/tree_based_stat_rmTS.updated.csv \
   FilterPara_rmTS/tree_based_stat_rmTS.csv

csvtk plot hist -f treeness --format pdf \
    -o FilterPara_rmTS/04_treeness_rmTS.pdf \
    < FilterPara_rmTS/tree_based_stat_rmTS.csv


# 3.2.5 Degree of violation of a molecular clock (DVMC)

TMP_DIR="FilterPara_rmTS/tmp5.0.4_dvmc"
OUT_CSV="${TMP_DIR}/dvmc.csv"

mkdir -p "$TMP_DIR"
echo "Locus,DVMC" > "$OUT_CSV"

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"

    if [[ -s "$TREE" ]]; then
        DVMC=$(phykit degree_of_violation_of_a_molecular_clock "$TREE" 2>/dev/null)
        echo "${LOCI},${DVMC}" >> "$OUT_CSV"
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done

csvtk join -H -f Locus \
    FilterPara_rmTS/tree_based_stat_rmTS.csv \
    "$OUT_CSV" \
    > FilterPara_rmTS/tree_based_stat_rmTS.updated.csv

mv FilterPara_rmTS/tree_based_stat_rmTS.updated.csv \
   FilterPara_rmTS/tree_based_stat_rmTS.csv

csvtk plot hist -f DVMC --format pdf \
    -o FilterPara_rmTS/05_DVMC_rmTS.pdf \
    < FilterPara_rmTS/tree_based_stat_rmTS.csv


# 3.2.6 Further inspection of loci with high treeness or DVMC

# Loci were flagged using treeness >= 0.3490 or DVMC >= 0.55.
awk -F',' 'NR > 1 && ($3 >= 0.3490 || $4 >= 0.55) {print $1}' \
    FilterPara_rmTS/tree_based_stat_rmTS.csv \
    > FilterPara_rmTS/high_trnsDVMC_loci.txt

echo "Locus,SeqID" > treeshrink_trnsDVMC.csv

MODE="per-gene"
QUANTILES="0.5"
TREE_DIR="genetrees_rmTS"
LOCI_LIST=$(cat FilterPara_rmTS/high_trnsDVMC_loci.txt)

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"
    TMP_DIR="tmp_treeshrink_trnsDVMC/${LOCI}"

    if [[ -s "$TREE" ]]; then
        mkdir -p "$TMP_DIR"
        cp "$TREE" "$TMP_DIR/input.tree"

        (
            cd "$TMP_DIR" || exit 1
            run_treeshrink.py -t input.tree -m "$MODE" -q "$QUANTILES" -O output > treeshrink.log 2>&1
        )

        RMSEQ_FILE="${TMP_DIR}/input_treeshrink/output.txt"
        if [[ -f "$RMSEQ_FILE" ]]; then
            while IFS=$'\t' read -r -a SEQS; do
                for ID in "${SEQS[@]}"; do
                    if [[ -n "$ID" ]]; then
                        echo "${LOCI},${ID}" >> treeshrink_trnsDVMC.csv
                    fi
                done
            done < "$RMSEQ_FILE"
        else
            echo "[Warning] TreeShrink output not found for $LOCI" >&2
        fi
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done

# Concatenate flagged gene trees for manual inspection.
OUT_LIST="FilterPara_rmTS/high_trnsDVMC_loci.txt"
OUT_NEXUS="concatenated_labeled_trees.nex"

{
    echo "#NEXUS"
    echo
    echo "Begin trees;"

    while read -r LOCUS; do
        TREE_FILE="${TREE_DIR}/${LOCUS}.treefile"
        if [[ -s "$TREE_FILE" ]]; then
            TREE_STRING=$(tr -d '\n' < "$TREE_FILE" | sed 's/;*$//')
            echo "  tree ${LOCUS} = ${TREE_STRING};"
        else
            echo "[Warning] Missing tree for $LOCUS" >&2
        fi
    done < "$OUT_LIST"

    echo "End;"
} > "$OUT_NEXUS"

# TreeShrink suggestions were checked manually. The curated removal table
# should be saved as treeshrink_trnsDVMC_manual.csv before this step.
python3 "$APPLY_SCRIPT" \
    -l "$LOCI_LIST_FILE" \
    -r treeshrink_trnsDVMC_manual.csv \
    -i ../1-raw_loci_rmTS \
    -o ../1-raw_loci_final \
    -c sequence_number_final.log \
    -n loci_name_final.log
