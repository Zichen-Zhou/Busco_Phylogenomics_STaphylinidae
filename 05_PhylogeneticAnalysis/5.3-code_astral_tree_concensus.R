# 5.3 Consensus across ASTRAL topologies

library(ape)

tree_files <- list.files(pattern = "\\.nwk$")
trees <- lapply(tree_files, read.tree)

# Confirm that all trees contain the same tip labels.
tip_sets <- lapply(trees, function(tree) sort(tree$tip.label))

if (length(tip_sets) > 1) {
  for (i in 2:length(tip_sets)) {
    if (!identical(tip_sets[[1]], tip_sets[[i]])) {
      stop(paste("Tip labels differ in tree:", tree_files[i]))
    }
  }
}

# Majority-rule consensus.
majority_consensus <- consensus(trees, p = 0.5)
write.tree(
  majority_consensus,
  file = "ASTRAL_14trees_majority_consensus.nwk"
)

# Strict consensus.
strict_consensus <- consensus(trees, p = 1.0)
write.tree(
  strict_consensus,
  file = "ASTRAL_14trees_strict_consensus.nwk"
)
