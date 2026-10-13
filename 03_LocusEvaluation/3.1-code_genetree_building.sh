# 3.1 Gene-tree building

# The same gene-tree inference procedure was used at three stages:
# 1) raw locus alignments, for preliminary long-branch sequence detection;
# 2) alignments after preliminary long-branch removal, for treeness- and DVMC-based checking;
# 3) alignments after non-ortholog removal, for coalescent-based phylogenetic analyses.
#
# The code below shows the SLURM-array implementation. Update the relative
# input and output directories for the corresponding gene-tree inference round.

# Example SLURM resources used:
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=1G
# #SBATCH --array=1-2124%128

LOCI_LIST="2-loci_filter/loci_name_filter.log"
ALIGN_DIR="4-loci_trim_magus/faa"
TREE_OUT_DIR="5-treebased/gene_tree"
LOCUS_WORK_DIR="5-treebased/gene_tree_work"
THREADS=1

mkdir -p "$TREE_OUT_DIR" "$LOCUS_WORK_DIR"

LOCI=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LOCI_LIST")
ALIGN="${ALIGN_DIR}/${LOCI}.aa.fas"
OUTDIR="${LOCUS_WORK_DIR}/${LOCI}"

if [[ ! -s "$ALIGN" ]]; then
    echo "[Warning] Missing alignment: $ALIGN" >&2
else
    mkdir -p "$OUTDIR"
    cp "$ALIGN" "$OUTDIR/"

    (
        cd "$OUTDIR" || exit 1
        iqtree2 -s "${LOCI}.aa.fas" -mset EX_EHO -bb 1000 -T "$THREADS" > "${LOCI}.log" 2>&1
    )

    TREEFILE="${OUTDIR}/${LOCI}.aa.fas.treefile"
    if [[ -s "$TREEFILE" ]]; then
        cp "$TREEFILE" "${TREE_OUT_DIR}/${LOCI}.treefile"
        echo "[OK] Tree completed for $LOCI"
    else
        echo "[Warning] Treefile missing for $LOCI" >&2
    fi
fi
