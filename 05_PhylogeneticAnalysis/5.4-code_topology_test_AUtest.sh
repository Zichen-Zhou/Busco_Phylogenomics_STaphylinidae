# 5.4 Topology tests using the approximately unbiased (AU) test
# FLT_taxa50 is shown as an example.
# Internal labels taxa50/75/90 correspond to matrix50/75/90 in the final manuscript.


# 5.4.1 Compare topologies inferred using different methods

cat *.treefile > trees4AU_FLT_taxa50.treels

iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m LG+C40+FO+R4 \
    -z trees4AU_FLT_taxa50.treels \
    -n 0 \
    -zb 2000 \
    -au \
    -pre au_FLT_taxa50 \
    -fs FcC_supermatrix_FLT_taxa50.fas.sitefreq \
    -nt AUTO


# 5.4.2 Infer trees under alternative backbone constraints

# Four-lineage hypothesis.
iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m LG+C60+FO+R \
    -fs FcC_supermatrix_FLT_taxa50.fas.sitefreq \
    -B 1000 \
    -alrt 1000 \
    -T 16 \
    -g 4group_backbone_w.og.nwk

# SSE-S hypothesis.
iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m LG+C60+FO+R \
    -fs FcC_supermatrix_FLT_taxa50.fas.sitefreq \
    -B 1000 \
    -alrt 1000 \
    -T 16 \
    -g SSE-S_backbone_w.og.nwk


# 5.4.3 AU test including hypothesis-constrained topologies

iqtree2 \
    -s FcC_supermatrix_FLT_taxa50.fas \
    -m LG+C60+FO+R \
    -z trees4AU_FLT_taxa50.treels \
    -n 0 \
    -zb 2000 \
    -au \
    -pre au_FLT_taxa50_PMSF \
    -fs FcC_supermatrix_FLT_taxa50.fas.sitefreq \
    -nt AUTO
