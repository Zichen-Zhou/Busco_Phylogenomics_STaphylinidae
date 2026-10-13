# 3.4 Locus evaluation
# Run each section from the project working directory unless otherwise noted.


# 3.4.1 TrimAl after non-ortholog removal

LOCI_NAME_FILTER=$(cat 2-loci_filter/loci_name_filter.log)

mkdir -p 4-loci_trim_final/fna 4-loci_trim_final/faa

for LOCI_FILTER in $LOCI_NAME_FILTER; do
    trimal \
        -in "3-faa_aln_final/${LOCI_FILTER}.faa" \
        -out "4-loci_trim_final/faa/${LOCI_FILTER}.aa.fas" \
        -automated1

    trimal \
        -in "3-faa_aln_final/${LOCI_FILTER}.faa" \
        -out "4-loci_trim_final/fna/${LOCI_FILTER}.nuc.fas" \
        -automated1 \
        -backtrans "1-raw_loci_final/fna/${LOCI_FILTER}.fna"
done


# 3.4.2 Tree-based evaluation

# 3.4.2.0 Gene-tree building
## See 3.1-code_genetree_building.sh.


# 3.4.2.1 Average branch support (ABS)

LOCI_LIST=$(cat 2-loci_filter/loci_name_filter.log)
TREE_DIR="6-dat.eva_treebased/genetrees_final.clean"
TREE_STAT="6-dat.eva_treebased/tree_based_stat.csv"

mkdir -p 6-dat.eva_treebased
echo "Locus,ABS" > "$TREE_STAT"

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"

    if [[ -s "$TREE" ]]; then
        ABS=$(sed 's/)/\n/g' "$TREE" \
            | cut -d ":" -f 1 \
            | sed '1d' \
            | sed 's/;//g' \
            | awk NF \
            | awk '{a+=$1} END {if (NR>0) print a/NR; else print "NA"}')
        echo "${LOCI},${ABS}" >> "$TREE_STAT"
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done

csvtk plot hist -f ABS --format pdf \
    -o 6-dat.eva_treebased/01_ABS.pdf \
    < "$TREE_STAT"


# 3.4.2.2 Treeness

TMP_DIR="6-dat.eva_treebased/tmp6.2.2_treeness"
TREENESS_CSV="${TMP_DIR}/treeness.csv"

mkdir -p "$TMP_DIR"
echo "Locus,treeness" > "$TREENESS_CSV"

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"

    if [[ -s "$TREE" ]]; then
        TREENESS=$(phykit treeness "$TREE" 2>/dev/null)
        echo "${LOCI},${TREENESS}" >> "$TREENESS_CSV"
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done

csvtk join -H -f Locus \
    "$TREE_STAT" "$TREENESS_CSV" \
    > 6-dat.eva_treebased/tree_based_stat.updated.csv

mv 6-dat.eva_treebased/tree_based_stat.updated.csv "$TREE_STAT"

csvtk plot hist -f treeness --format pdf \
    -o 6-dat.eva_treebased/02_treeness.pdf \
    < "$TREE_STAT"


# TOR requires the RCV values generated in section 3.4.3.2.
# Run this block after completing section 3.4.3.2.

csvtk join -H -f Locus \
    6-dat.eva_alnbased/aln_based_stat.csv \
    "$TREE_STAT" \
    > tmp_data.eva_joined.csv

# In the joined table used in this study, RCV was column 10 and treeness
# was column 15. TOR = treeness / RCV.
awk -F',' \
    'NR==1 {print $0",TOR"; next} {printf "%s,%.6f\n", $0, ($15/$10)}' \
    tmp_data.eva_joined.csv \
    > data.eva_stat.csv

rm tmp_data.eva_joined.csv

csvtk join -f Locus \
    "$TREE_STAT" \
    <(csvtk cut -f Locus,TOR data.eva_stat.csv) \
    > 6-dat.eva_treebased/tree_based_stat.updated.csv

mv data.eva_stat.csv \
   6-dat.eva_treebased/tmp6.2.2_treeness/data.eva_stat_TOR.csv

mv 6-dat.eva_treebased/tree_based_stat.updated.csv "$TREE_STAT"

csvtk plot hist -f TOR --format pdf \
    -o 6-dat.eva_treebased/02_TOR.pdf \
    < "$TREE_STAT"

csvtk filter -f "TOR<=1" "$TREE_STAT" \
    | csvtk plot hist -f TOR --bins 100 --format pdf \
        -o 6-dat.eva_treebased/02_TOR_adjusted.pdf


# 3.4.2.3 Degree of violation of a molecular clock (DVMC)

TMP_DIR="6-dat.eva_treebased/tmp6.2.3_DVMC"
DVMC_CSV="${TMP_DIR}/dvmc.csv"

mkdir -p "$TMP_DIR"
echo "Locus,DVMC" > "$DVMC_CSV"

for LOCI in $LOCI_LIST; do
    TREE="${TREE_DIR}/${LOCI}.treefile"

    if [[ -s "$TREE" ]]; then
        DVMC=$(phykit degree_of_violation_of_a_molecular_clock "$TREE" 2>/dev/null)
        echo "${LOCI},${DVMC}" >> "$DVMC_CSV"
    else
        echo "[Warning] Missing tree for $LOCI" >&2
    fi
done

csvtk join -H -f Locus \
    "$TREE_STAT" "$DVMC_CSV" \
    > 6-dat.eva_treebased/tree_based_stat.updated.csv

mv 6-dat.eva_treebased/tree_based_stat.updated.csv "$TREE_STAT"

csvtk plot hist -f DVMC --format pdf \
    -o 6-dat.eva_treebased/03_DVMC.pdf \
    < "$TREE_STAT"

csvtk filter -f "DVMC<=1" "$TREE_STAT" \
    | csvtk plot hist -f DVMC --bins 100 --format pdf \
        -o 6-dat.eva_treebased/03_DVMC_adjusted.pdf


# 3.4.3 Alignment-based evaluation

mkdir -p 6-dat.eva_alnbased
ALN_STAT="6-dat.eva_alnbased/aln_based_stat.csv"
LOCI_NAME_FILTER=$(cat 2-loci_filter/loci_name_filter.log)


# 3.4.3.1 Alignment length and parsimony-informative sites

echo "Locus,Parsimony_sites_nt,Length_nt,Percentage_nt" > "$ALN_STAT"

for LOCI in $LOCI_NAME_FILTER; do
    INPUT_FILE="4-loci_trim_final.clean/fna/${LOCI}.nuc.fas"

    if [[ -s "$INPUT_FILE" ]]; then
        RESULT=$(phykit pis "$INPUT_FILE" 2>/dev/null)
        PARS=$(echo "$RESULT" | awk '{print $1}')
        LEN=$(echo "$RESULT" | awk '{print $2}')
        PCT=$(echo "$RESULT" | awk '{print $3}')
        echo "$LOCI,$PARS,$LEN,$PCT" >> "$ALN_STAT"
    else
        echo "[Warning] Missing or empty input: $INPUT_FILE" >&2
    fi
done

csvtk plot hist -f Parsimony_sites_nt --format pdf \
    -o 6-dat.eva_alnbased/01_Parsimony_sites_nt.pdf \
    < "$ALN_STAT"

csvtk plot hist -f Length_nt --format pdf \
    -o 6-dat.eva_alnbased/01_Alignment_length_nt.pdf \
    < "$ALN_STAT"

csvtk plot hist -f Percentage_nt --format pdf \
    -o 6-dat.eva_alnbased/01_Percentage_of_parsimony_sites_nt.pdf \
    < "$ALN_STAT"

TMP_DIR="6-dat.eva_alnbased/tmp_6.1PisLenPer"
TMP_AA_STAT="${TMP_DIR}/tmp_aa_stat.csv"

mkdir -p "$TMP_DIR"
echo "Locus,Parsimony_sites_aa,Length_aa,Percentage_aa" > "$TMP_AA_STAT"

for LOCI in $LOCI_NAME_FILTER; do
    INPUT_FILE="4-loci_trim_final.clean/faa/${LOCI}.aa.fas"

    if [[ -s "$INPUT_FILE" ]]; then
        RESULT=$(phykit pis "$INPUT_FILE" 2>/dev/null)
        PARS=$(echo "$RESULT" | awk '{print $1}')
        LEN=$(echo "$RESULT" | awk '{print $2}')
        PCT=$(echo "$RESULT" | awk '{print $3}')
        echo "$LOCI,$PARS,$LEN,$PCT" >> "$TMP_AA_STAT"
    else
        echo "[Warning] Missing or empty input: $INPUT_FILE" >&2
    fi
done

csvtk join -H -f Locus \
    "$ALN_STAT" "$TMP_AA_STAT" \
    > 6-dat.eva_alnbased/aln_based_stat.updated.csv

mv 6-dat.eva_alnbased/aln_based_stat.updated.csv "$ALN_STAT"

csvtk plot hist -f Parsimony_sites_aa --format pdf \
    -o 6-dat.eva_alnbased/01_Parsimony_sites_aa.pdf \
    < "$ALN_STAT"

csvtk plot hist -f Length_aa --format pdf \
    -o 6-dat.eva_alnbased/01_Alignment_length_aa.pdf \
    < "$ALN_STAT"

csvtk plot hist -f Percentage_aa --format pdf \
    -o 6-dat.eva_alnbased/01_Percentage_of_parsimony_sites_aa.pdf \
    < "$ALN_STAT"


# 3.4.3.2 GC content and relative composition variability (RCV)

GC_DIR="6-dat.eva_alnbased/tmp_6.2GC"
RCV_DIR="6-dat.eva_alnbased/tmp_6.2RCV"

mkdir -p "$GC_DIR" "$RCV_DIR"
echo "Locus,GC_content" > "${GC_DIR}/gc_stat.csv"
echo "Locus,RCV" > "${RCV_DIR}/rcv_stat.csv"

for LOCI in $LOCI_NAME_FILTER; do
    INPUT_FILE="4-loci_trim_final.clean/fna/${LOCI}.nuc.fas"

    if [[ -s "$INPUT_FILE" ]]; then
        GC=$(phykit gc_content "$INPUT_FILE" | awk '{print $1}')
        echo "${LOCI},${GC}" >> "${GC_DIR}/gc_stat.csv"
    else
        echo "[Warning] Missing or empty input: $INPUT_FILE" >&2
    fi
done

for LOCI in $LOCI_NAME_FILTER; do
    INPUT_FILE="4-loci_trim_final.clean/faa/${LOCI}.aa.fas"

    if [[ -s "$INPUT_FILE" ]]; then
        RCV=$(phykit relative_composition_variability "$INPUT_FILE" | awk '{print $1}')
        echo "${LOCI},${RCV}" >> "${RCV_DIR}/rcv_stat.csv"
    else
        echo "[Warning] Missing or empty input: $INPUT_FILE" >&2
    fi
done

csvtk join -H -f Locus \
    "$ALN_STAT" \
    <(csvtk join -H -f Locus \
        "${GC_DIR}/gc_stat.csv" \
        "${RCV_DIR}/rcv_stat.csv") \
    > 6-dat.eva_alnbased/aln_based_stat.updated.csv

mv 6-dat.eva_alnbased/aln_based_stat.updated.csv "$ALN_STAT"

csvtk plot hist -f GC_content --format pdf \
    -o 6-dat.eva_alnbased/02_GC_content.pdf \
    < "$ALN_STAT"

csvtk plot hist -f RCV --format pdf \
    -o 6-dat.eva_alnbased/02_RCV.pdf \
    < "$ALN_STAT"


# 3.4.3.3 Symmetry test for stationarity, reversibility and homogeneity

SRH_DIR="6-dat.eva_alnbased/tmp_6.3SRH"
mkdir -p "$SRH_DIR"

(
    cd "$SRH_DIR" || exit 1

    rm -rf faa
    cp -r ../../4-loci_trim_final.clean/faa ./
    ls faa | sed 's|^|faa/|' > all.alignments

    phykit cc -a all.alignments -p all
    iqtree2 -s all.fa -p all.partition --symtest-only -T 4

    grep -v "^#" all.partition.symtest.csv \
        | awk -F',' 'NR>1 {
            fname = $1
            sub(/^faa\//, "", fname)
            sub(/\.aa\.fas$/, "", fname)
            print fname "," $4
        }' \
        > symtest.pvalues.txt

    sed -i '1iLocus,symtest_pval' symtest.pvalues.txt
)

csvtk join -H -f Locus \
    "$ALN_STAT" "${SRH_DIR}/symtest.pvalues.txt" \
    > 6-dat.eva_alnbased/aln_based_stat.updated.csv

mv 6-dat.eva_alnbased/aln_based_stat.updated.csv "$ALN_STAT"

csvtk plot hist -f symtest_pval --format pdf \
    -o 6-dat.eva_alnbased/03_symtest_pval.pdf \
    < "$ALN_STAT"


# 3.4.3.4 Average pairwise sequence identity

RATE_DIR="6-dat.eva_alnbased/tmp_6.4EvoRate"
mkdir -p "$RATE_DIR"

echo "Locus,evolutionary_rate_nt" > "${RATE_DIR}/evo_rate_nt.txt"

for LOCI in $LOCI_NAME_FILTER; do
    INPUT_FILE="4-loci_trim_final.clean/fna/${LOCI}.nuc.fas"

    if [[ -s "$INPUT_FILE" ]]; then
        RATE=$(phykit pairwise_identity "$INPUT_FILE" 2>/dev/null \
            | grep "^mean:" \
            | awk '{print $2}')
        echo "${LOCI},${RATE}" >> "${RATE_DIR}/evo_rate_nt.txt"
    else
        echo "[Warning] Missing or empty input: $INPUT_FILE" >&2
    fi
done

csvtk join -H -f Locus \
    "$ALN_STAT" "${RATE_DIR}/evo_rate_nt.txt" \
    > 6-dat.eva_alnbased/aln_based_stat.updated.csv

mv 6-dat.eva_alnbased/aln_based_stat.updated.csv "$ALN_STAT"

csvtk plot hist -f evolutionary_rate_nt --format pdf \
    -o 6-dat.eva_alnbased/04_EvoRate_nt_hist.pdf \
    < "${RATE_DIR}/evo_rate_nt.txt"

echo "Locus,evolutionary_rate_aa" > "${RATE_DIR}/evo_rate_aa.txt"

for LOCI in $LOCI_NAME_FILTER; do
    INPUT_FILE="4-loci_trim_final.clean/faa/${LOCI}.aa.fas"

    if [[ -s "$INPUT_FILE" ]]; then
        RATE=$(phykit pairwise_identity "$INPUT_FILE" 2>/dev/null \
            | grep "^mean:" \
            | awk '{print $2}')
        echo "${LOCI},${RATE}" >> "${RATE_DIR}/evo_rate_aa.txt"
    else
        echo "[Warning] Missing or empty input: $INPUT_FILE" >&2
    fi
done

csvtk join -H -f Locus \
    "$ALN_STAT" "${RATE_DIR}/evo_rate_aa.txt" \
    > 6-dat.eva_alnbased/aln_based_stat.updated.csv

mv 6-dat.eva_alnbased/aln_based_stat.updated.csv "$ALN_STAT"

csvtk plot hist -f evolutionary_rate_aa --format pdf \
    -o 6-dat.eva_alnbased/04_EvoRate_aa_hist.pdf \
    < "${RATE_DIR}/evo_rate_aa.txt"


# 3.4.3.5 Likelihood mapping

# The following block was run as a SLURM array.
# Example SLURM resources used:
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=250M
# #SBATCH --array=1-2124%128

ALIGN_DIR="4-loci_trim_final.clean/faa"
LMAP_DIR="6-dat.eva_alnbased/tmp_6.5lmap"

mkdir -p "$LMAP_DIR"

LOCI=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 2-loci_filter/loci_name_filter.log)
ALIGN="${ALIGN_DIR}/${LOCI}.aa.fas"
OUT_PREFIX="${LMAP_DIR}/${LOCI}"

if [[ ! -s "$ALIGN" ]]; then
    echo "[Warning] Missing alignment: $ALIGN" >&2
    echo "${LOCI},NA" > "${LMAP_DIR}/${LOCI}.lmap.tmp"
else
    iqtree2 \
        -s "$ALIGN" \
        -m Q.insect \
        -lmap ALL \
        -n 0 \
        -pre "$OUT_PREFIX" \
        -quiet

    IQTREE_FILE="${OUT_PREFIX}.iqtree"

    if [[ -f "$IQTREE_FILE" ]]; then
        VALUE=$(grep "Number of fully resolved  quartets" "$IQTREE_FILE" \
            | sed -n 's/.*(=\([0-9.]\+%\)).*/\1/p')
        [[ -z "$VALUE" ]] && VALUE="NA"
    else
        VALUE="NA"
    fi

    echo "${LOCI},${VALUE}" > "${LMAP_DIR}/${LOCI}.lmap.tmp"
fi


# Summarize likelihood-mapping results after all array jobs are complete.

LMAP_STAT="${LMAP_DIR}/lmap_stat.csv"
echo "Locus,FullyResolvedPct" > "$LMAP_STAT"
cat "${LMAP_DIR}"/*.lmap.tmp >> "$LMAP_STAT"

csvtk join -H -f Locus \
    "$ALN_STAT" "$LMAP_STAT" \
    > 6-dat.eva_alnbased/aln_based_stat.updated.csv

mv 6-dat.eva_alnbased/aln_based_stat.updated.csv "$ALN_STAT"

# Convert percentages to decimal values while retaining missing values as NA.
awk -F',' '
    BEGIN {OFS=","}
    NR==1 {print; next}
    {
        if ($NF == "" || $NF == "NA") {
            print
            next
        }
        gsub(/%/, "", $NF)
        $NF = $NF / 100
        print
    }
' "$ALN_STAT" \
    > 6-dat.eva_alnbased/aln_based_stat.lmap_decimal.csv

# Exclude missing values only from the histogram.
awk -F',' 'NR==1 || ($NF != "" && $NF != "NA")' \
    6-dat.eva_alnbased/aln_based_stat.lmap_decimal.csv \
    | csvtk plot hist -f FullyResolvedPct --format pdf \
        -o 6-dat.eva_alnbased/05_Likelihood_mapping.pdf
