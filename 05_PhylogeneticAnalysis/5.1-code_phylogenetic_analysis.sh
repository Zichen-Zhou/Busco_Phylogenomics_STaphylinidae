# 5.1 Phylogenetic analyses
# FLT_taxa50 is shown as an example.
# Internal labels taxa50/75/90 correspond to matrix50/75/90 in the final manuscript.


# 5.1.1 LG

iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -p FcC_supermatrix_partition_FLT_taxa50.txt \
    -m MFP \
    --mset LG \
    --msub nuclear \
    --rclusterf 10 \
    -B 1000 \
    --alrt 1000 \
    -T AUTO


# 5.1.2 GHOST

iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m LG+FO+H4 \
    -B 1000 \
    --alrt 1000 \
    -T 3 \
    -safe

python3 5.1-script_extract_main_GHOSTtree.py \
    FcC_supermatrix_FLT_taxa50.fas.treefile \
    GHOST_FLT_taxa50.nwk


# 5.1.3 EX_EHO

iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m EX_EHO+FO+R \
    -B 1000 \
    --alrt 1000 \
    -T 16 \
    --safe


# 5.1.4 PMSF

# Step 1: estimate site-specific amino-acid profiles.
iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m LG+C40+FO+R4 \
    -ft LG_FLT_taxa50.txt.treefile \
    -n 0 \
    -T 8

# Step 2: infer the PMSF tree using the estimated site-frequency file.
iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m LG+C60+FO+R \
    -fs FcC_supermatrix_FLT_taxa50.fas.sitefreq \
    -B 1000 \
    -alrt 1000 \
    -T 16


# 5.1.5 Dayhoff6 recoding

perl pgrecodeseq.pl \
    --type=ANY \
    "ARNDCQEGHILKMFPSTWYVX-01223220144145000554?" \
    FcC_supermatrix_FLT_taxa50.fas \
    Dayhoff6_FLT_taxa50.fas

iqtree2 \
    -s Dayhoff6_FLT_taxa50.fas \
    -m GTR+R \
    -B 1000 \
    -alrt 1000 \
    -T 8


# 5.1.6 ASTRAL

TREE_DIR="6-dat.eva_treebased_len70/genetrees_final.clean_len70"
LIST="loci_list_FLT_taxa50.txt"
OUT="genetrees_FLT_taxa50.tre"

: > "$OUT"

sed 's/\r$//' "$LIST" | while IFS= read -r LOCUS; do
    [[ -z "$LOCUS" ]] && continue
    cat "${TREE_DIR}/${LOCUS}.treefile" >> "$OUT"
done

astral \
    -i genetrees_FLT_taxa50.tre \
    -o ASTRAL_FLT_taxa50.LPP.treefile \
    2> ASTRAL_FLT_taxa50.LPP.log


# Calculate gCF and sCF.
# The same IQ-TREE concordance-factor procedure was used for concatenation-based trees.

ALIGN_DIR="4-loci_trim_final.clean_len70/faa"
LIST="loci_list_FLT_taxa50.txt"
ALN_SUBDIR="aln_FLT_taxa50_subset"

mkdir -p "$ALN_SUBDIR"

sed 's/\r$//' "$LIST" | while IFS= read -r LOCUS; do
    [[ -z "$LOCUS" ]] && continue
    ln -sfn \
        "${ALIGN_DIR}/${LOCUS}_len70.aa.fas" \
        "${ALN_SUBDIR}/${LOCUS}_len70.aa.fas"
done

iqtree2 \
    -t ASTRAL_FLT_taxa50.LPP.treefile \
    --gcf genetrees_FLT_taxa50.tre \
    -p "$ALN_SUBDIR" \
    --scf 100 \
    -pre ASTRAL_concord_FLT_taxa50 \
    -T 8
