#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

# ============================================================
# AU candidate-tree topology deduplication
#
# Functions:
#   1. Read multiple candidate trees from a .treels file
#   2. Remove branch lengths
#   3. Remove internal node labels/support values
#   4. Convert trees to unrooted topology
#   5. Deduplicate identical unrooted topologies
#   6. Write:
#        - topology-only unique trees for IQ-TREE AU test
#        - original-tree -> topology mapping
#        - unique-topology summary
#
# Usage:
#
# Rscript dedup_AU_trees.R \
#   -i input.treels \
#   -m LG,GHOST,EXEHO,PMSF,Dayhoff6,ASTRAL \
#   -o output_directory \
#   -p prefix
#
# -m is optional. If omitted, trees are called Tree1, Tree2, ...
# ============================================================


# -------------------------
# Help
# -------------------------

show_help <- function() {
  
  cat("
Usage:

  Rscript dedup_AU_trees.R -i INPUT [options]

Required:

  -i, --input FILE
      Input file containing multiple Newick trees.

Optional:

  -m, --methods LIST
      Comma-separated method/model names in exactly the same order
      as the trees in the input file.

      Example:
      -m LG,GHOST,EXEHO,PMSF,Dayhoff6,ASTRAL

      If omitted, methods will be named:
      Tree1, Tree2, Tree3, ...

  -o, --outdir DIR
      Output directory.
      Default: current directory

  -p, --prefix STRING
      Prefix for output files.
      Default: basename of input file without .treels

  -h, --help
      Show this help message.

Example:

  Rscript dedup_AU_trees.R \\
      -i trees4AU_FLT_taxa50.treels \\
      -m LG,GHOST,EXEHO,PMSF,Dayhoff6,ASTRAL \\
      -o AU_dedup \\
      -p FLT_taxa50

Outputs:

  PREFIX.topology_unique.treels
  PREFIX.topology_mapping.tsv
  PREFIX.unique_topologies.tsv

")
}


# -------------------------
# Parse command line
# -------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 ||
    any(args %in% c("-h", "--help"))) {
  
  show_help()
  quit(status = 0)
}

input_file <- NULL
method_string <- NULL
outdir <- "."
prefix <- NULL

i <- 1

while (i <= length(args)) {
  
  arg <- args[i]
  
  if (arg %in% c("-i", "--input")) {
    
    i <- i + 1
    
    if (i > length(args)) {
      stop("Missing value after ", arg)
    }
    
    input_file <- args[i]
    
  } else if (arg %in% c("-m", "--methods")) {
    
    i <- i + 1
    
    if (i > length(args)) {
      stop("Missing value after ", arg)
    }
    
    method_string <- args[i]
    
  } else if (arg %in% c("-o", "--outdir")) {
    
    i <- i + 1
    
    if (i > length(args)) {
      stop("Missing value after ", arg)
    }
    
    outdir <- args[i]
    
  } else if (arg %in% c("-p", "--prefix", "-prefix")) {
    
    i <- i + 1
    
    if (i > length(args)) {
      stop("Missing value after ", arg)
    }
    
    prefix <- args[i]
    
  } else {
    
    stop("Unknown argument: ", arg)
  }
  
  i <- i + 1
}


# -------------------------
# Validate input
# -------------------------

if (is.null(input_file)) {
  stop("Input file is required. Use -i FILE")
}

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

if (!dir.exists(outdir)) {
  
  dir.create(
    outdir,
    recursive = TRUE
  )
}

if (is.null(prefix)) {
  
  prefix <- basename(input_file)
  
  prefix <- sub(
    "\\.treels$",
    "",
    prefix,
    ignore.case = TRUE
  )
}


# -------------------------
# Read trees
# -------------------------

trees <- read.tree(input_file)

if (inherits(trees, "phylo")) {
  
  trees <- structure(
    list(trees),
    class = "multiPhylo"
  )
}

n_trees <- length(trees)

cat("\nInput file: ", input_file, "\n", sep = "")
cat("Number of input trees: ", n_trees, "\n", sep = "")


# -------------------------
# Method names
# -------------------------

if (is.null(method_string)) {
  
  methods <- paste0(
    "Tree",
    seq_len(n_trees)
  )
  
} else {
  
  methods <- strsplit(
    method_string,
    ",",
    fixed = TRUE
  )[[1]]
  
  methods <- trimws(methods)
  
  if (length(methods) != n_trees) {
    
    stop(
      "\nNumber of methods does not match number of trees.\n",
      "Trees:   ", n_trees, "\n",
      "Methods: ", length(methods), "\n"
    )
  }
  
  if (any(methods == "")) {
    stop("Empty method name detected in -m argument.")
  }
  
  if (anyDuplicated(methods)) {
    warning(
      "Duplicated method names were supplied. ",
      "This is allowed, but check that the order is correct."
    )
  }
}


# -------------------------
# Check tip labels
# -------------------------

reference_tips <- sort(
  trees[[1]]$tip.label
)

for (j in seq_len(n_trees)) {
  
  current_tips <- sort(
    trees[[j]]$tip.label
  )
  
  if (!identical(
    reference_tips,
    current_tips
  )) {
    
    missing_tips <- setdiff(
      reference_tips,
      current_tips
    )
    
    extra_tips <- setdiff(
      current_tips,
      reference_tips
    )
    
    stop(
      "\nTree ", j,
      " (", methods[j], ") ",
      "does not contain the same tip set.\n",
      "Missing: ",
      paste(missing_tips, collapse = ", "),
      "\nExtra: ",
      paste(extra_tips, collapse = ", "),
      "\n"
    )
  }
}

cat(
  "All trees contain the same ",
  length(reference_tips),
  " tips.\n",
  sep = ""
)


# -------------------------
# Convert to pure
# unrooted topologies
# -------------------------

topology_trees <- lapply(
  trees,
  function(tr) {
    
    # Remove branch lengths
    tr$edge.length <- NULL
    
    # Remove internal support/node labels
    tr$node.label <- NULL
    
    # Convert to unrooted topology
    if (is.rooted(tr)) {
      tr <- unroot(tr)
    }
    
    # Standardize edge ordering
    tr <- reorder.phylo(
      tr,
      order = "cladewise"
    )
    
    tr
  }
)

class(topology_trees) <- "multiPhylo"


# -------------------------
# Pairwise topology distance
# -------------------------

rf <- as.matrix(
  dist.topo(
    topology_trees,
    method = "PH85"
  )
)

rownames(rf) <- paste0(
  "Tree",
  seq_len(n_trees)
)

colnames(rf) <- rownames(rf)


# -------------------------
# Deduplicate
# -------------------------

topology_id <- rep(
  NA_character_,
  n_trees
)

au_tree_index <- rep(
  NA_integer_,
  n_trees
)

representatives <- integer(0)

for (j in seq_len(n_trees)) {
  
  if (length(representatives) == 0) {
    
    representatives <- j
    
    topology_id[j] <- "AU_T01"
    au_tree_index[j] <- 1
    
    next
  }
  
  identical_reps <- representatives[
    rf[j, representatives] == 0
  ]
  
  if (length(identical_reps) > 0) {
    
    representative <- identical_reps[1]
    
    topology_number <- match(
      representative,
      representatives
    )
    
    topology_id[j] <- sprintf(
      "AU_T%02d",
      topology_number
    )
    
    au_tree_index[j] <- topology_number
    
  } else {
    
    representatives <- c(
      representatives,
      j
    )
    
    topology_number <- length(
      representatives
    )
    
    topology_id[j] <- sprintf(
      "AU_T%02d",
      topology_number
    )
    
    au_tree_index[j] <- topology_number
  }
}


# -------------------------
# Mapping table
# -------------------------

mapping <- data.frame(
  
  original_tree_index = seq_len(n_trees),
  
  original_tree = paste0(
    "Tree",
    seq_len(n_trees)
  ),
  
  method = methods,
  
  AU_tree_index = au_tree_index,
  
  topology_id = topology_id,
  
  representative = seq_len(n_trees) %in%
    representatives,
  
  stringsAsFactors = FALSE
)


# Add representative details

mapping$representative_original_tree <- NA_character_
mapping$representative_method <- NA_character_

for (j in seq_len(n_trees)) {
  
  rep_index <- representatives[
    mapping$AU_tree_index[j]
  ]
  
  mapping$representative_original_tree[j] <-
    paste0("Tree", rep_index)
  
  mapping$representative_method[j] <-
    methods[rep_index]
}


# -------------------------
# Unique topology summary
# -------------------------

unique_summary <- do.call(
  
  rbind,
  
  lapply(
    seq_along(representatives),
    
    function(k) {
      
      members <- which(
        au_tree_index == k
      )
      
      rep_index <- representatives[k]
      
      data.frame(
        
        AU_tree_index = k,
        
        topology_id = sprintf(
          "AU_T%02d",
          k
        ),
        
        representative_original_tree =
          paste0(
            "Tree",
            rep_index
          ),
        
        representative_method =
          methods[rep_index],
        
        original_trees =
          paste(
            paste0(
              "Tree",
              members
            ),
            collapse = ";"
          ),
        
        methods =
          paste(
            methods[members],
            collapse = ";"
          ),
        
        n_original_trees =
          length(members),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


# -------------------------
# Unique tree set
# -------------------------

unique_trees <- topology_trees[
  representatives
]


# -------------------------
# Output filenames
# -------------------------

tree_out <- file.path(
  outdir,
  paste0(
    prefix,
    ".topology_unique.treels"
  )
)

mapping_out <- file.path(
  outdir,
  paste0(
    prefix,
    ".topology_mapping.tsv"
  )
)

summary_out <- file.path(
  outdir,
  paste0(
    prefix,
    ".unique_topologies.tsv"
  )
)


# -------------------------
# Write outputs
# -------------------------

write.tree(
  unique_trees,
  file = tree_out
)

write.table(
  mapping,
  file = mapping_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  unique_summary,
  file = summary_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# -------------------------
# Report
# -------------------------

cat("\nOriginal tree mapping:\n\n")

print(
  mapping,
  row.names = FALSE
)

cat("\nUnique topology summary:\n\n")

print(
  unique_summary,
  row.names = FALSE
)

cat(
  "\n",
  n_trees,
  " input trees -> ",
  length(representatives),
  " unique unrooted topologies\n",
  sep = ""
)

cat("\nOutput files:\n")

cat(
  "  ",
  tree_out,
  "\n",
  sep = ""
)

cat(
  "  ",
  mapping_out,
  "\n",
  sep = ""
)

cat(
  "  ",
  summary_out,
  "\n",
  sep = ""
)

