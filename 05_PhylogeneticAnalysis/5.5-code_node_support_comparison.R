# 5.5 Node-support and concordance-factor comparison
# Concatenation-based trees are shown here; the same framework was also applied to ASTRAL trees.

# 5.5.1 Distribution of node-support and concordance-factor values
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Input and output
tree_files <- list.files(
  "0_treefiles_CF_rl_rt/",
  pattern = "\\.nwk$|\\.treefile$|\\.tre$|\\.cf\\.tree$",
  full.names = TRUE
)

dir.create("1_violin_CF", showWarnings = FALSE)

# Result containers
summary_list <- list()
plot_list    <- list()

# Process each tree
for (file in tree_files) {
  
  tr <- read.tree(file)
  
  node_labels <- tr$node.label
  node_labels <- node_labels[!is.na(node_labels) & node_labels != ""]
  
  split_labels <- strsplit(node_labels, "/")
  
  df <- do.call(rbind, lapply(split_labels, function(x) {
    if (length(x) == 4) {
      data.frame(
        UFboot  = as.numeric(x[1]),
        SH_aLRT = as.numeric(x[2]),
        gCF     = as.numeric(x[3]),
        sCF     = as.numeric(x[4])
      )
    } else NULL
  }))
  
  if (is.null(df) || nrow(df) == 0) next
  
  ## clip to [0,100]
  df$UFboot  <- pmin(pmax(df$UFboot,  0), 100)
  df$SH_aLRT <- pmin(pmax(df$SH_aLRT, 0), 100)
  df$gCF     <- pmin(pmax(df$gCF,     0), 100)
  df$sCF     <- pmin(pmax(df$sCF,     0), 100)
  
  ## summary statistics
  mean_vals <- colMeans(df, na.rm = TRUE)
  
  summary_list[[file]] <- data.frame(
    Tree = basename(file),
    Mean_UFboot = mean_vals["UFboot"],
    Mean_SHaLRT = mean_vals["SH_aLRT"],
    Mean_gCF    = mean_vals["gCF"],
    Mean_sCF    = mean_vals["sCF"],
    Max_UFboot  = max(df$UFboot, na.rm = TRUE),
    Max_SHaLRT  = max(df$SH_aLRT, na.rm = TRUE),
    Max_gCF     = max(df$gCF,    na.rm = TRUE),
    Max_sCF     = max(df$sCF,    na.rm = TRUE),
    Min_UFboot  = min(df$UFboot, na.rm = TRUE),
    Min_SHaLRT  = min(df$SH_aLRT, na.rm = TRUE),
    Min_gCF     = min(df$gCF,    na.rm = TRUE),
    Min_sCF     = min(df$sCF,    na.rm = TRUE),
    N_sCF_lt_33 = sum(df$sCF < 33.33, na.rm = TRUE)
  )
  
  df_long <- df %>%
    pivot_longer(
      cols = c("UFboot", "SH_aLRT", "gCF", "sCF"),
      names_to  = "Metric",
      values_to = "Value"
    )
  
  mean_df <- data.frame(
    Metric = c("UFboot", "SH_aLRT", "gCF", "sCF"),
    Mean   = c(
      mean_vals["UFboot"],
      mean_vals["SH_aLRT"],
      mean_vals["gCF"],
      mean_vals["sCF"]
    )
  )
  
  mean_df$Metric <- factor(
    mean_df$Metric,
    levels = c("UFboot", "SH_aLRT", "gCF", "sCF")
  )
  mean_df$x_center <- as.numeric(mean_df$Metric)
  
  ## VIOLIN PLOT
  p <- ggplot(
    df_long,
    aes(
      x    = factor(Metric, levels = c("UFboot", "SH_aLRT", "gCF", "sCF")),
      y    = Value,
      fill = Metric
    )
  ) +
    geom_violin(trim = TRUE, alpha = 0.7, color = "black") +
    
    # Mean line.
    geom_segment(
      data = mean_df,
      aes(
        x    = x_center - 0.32,
        xend = x_center + 0.32,
        y    = Mean,
        yend = Mean
      ),
      inherit.aes = FALSE,
      color = "red",
      size  = 0.6
    ) +
    
    # Mean label.
    geom_text(
      data = mean_df,
      aes(
        x     = x_center + 0.42,
        y     = Mean,
        label = sprintf("%.1f", Mean)
      ),
      inherit.aes = FALSE,
      color = "red",
      size  = 3,
      hjust = 0
    ) +
    
    scale_fill_manual(
      values = c(
        "UFboot"  = "firebrick",
        "SH_aLRT" = "gold3",
        "gCF"     = "#8499ce",
        "sCF"     = "#574774"
      ),
      breaks = c("UFboot", "SH_aLRT", "gCF", "sCF"),
      labels = c("UFboot", "SH-aLRT", "gCF", "sCF")
    ) +
    
    # Horizontal grid lines at 0/25/50/75/100.
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      minor_breaks = NULL
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      plot.margin = margin(10, 40, 10, 10)
    ) +
    
    labs(
      title = basename(file),
      x     = "Node support type",
      y     = "Support value"
    )
  
  outname <- paste0(
    "1_violin_CF/",
    tools::file_path_sans_ext(basename(file)),
    "_violin.png"
  )
  
  plot_list[[basename(file)]] <- p
  ggsave(outname, plot = p, width = 7, height = 5, dpi = 300)
}

# Summary table
summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, "1_violin_CF/node_support_summary_CF.csv", row.names = FALSE)

message("Violin plots with sCF generated successfully.")

# Combined panel
dataset_names <- c(
  "ALL_matrix50","ALL_matrix50_len70","ALL_matrix75",
  "ALL_matrix75_len70","ALL_matrix90","ALL_matrix90_len70",
  "SRH_matrix50", "SRH_matrix50_len70","SRH_matrix75",
  "SRH_matrix75_len70","SRH_matrix90","SRH_matrix90_len70",
  "FLT_matrix50", "FLT_matrix75"
)

method_levels <- c("LG", "GHOST", "EXEHO", "PMSF", "Dayhoff6")

plot_info <- data.frame(
  file = names(plot_list),
  base = tools::file_path_sans_ext(names(plot_list)),
  stringsAsFactors = FALSE
)

plot_info$Method  <- sub("_.*$", "", plot_info$base)
plot_info$Dataset <- NA_character_

for (ds in dataset_names) {
  plot_info$Dataset[grepl(ds, plot_info$base, fixed = TRUE)] <- ds
}

plot_info <- subset(plot_info,
                    !is.na(Dataset) & Method %in% method_levels)

plot_info$Method  <- factor(plot_info$Method,  levels = method_levels)
plot_info$Dataset <- factor(plot_info$Dataset, levels = dataset_names)

row_panels <- list()
for (ds in dataset_names) {
  sub_df <- subset(plot_info, Dataset == ds)
  row_plots <- vector("list", length(method_levels))
  names(row_plots) <- method_levels
  
  for (m in method_levels) {
    this <- sub_df[sub_df$Method == m, ]
    row_plots[[m]] <- if (nrow(this) == 1) plot_list[[this$file]] else ggplot() + theme_void()
  }
  row_panels[[ds]] <- wrap_plots(row_plots, ncol = length(method_levels))
}

final_panel <- wrap_plots(row_panels, ncol = 1)

ggsave(
  filename = "1_violin_CF/violin_panel_all_datasets_CF.png",
  plot     = final_panel,
  width    = 30,
  height   = 4 * length(dataset_names),
  dpi      = 300,
  limitsize = FALSE
)


# 5.5.2 gCF by hierarchical level
library(ape)
library(phangorn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

dir.create("gCF.vs.node", showWarnings = FALSE)

tree_files <- list.files(
  "0_treefiles_CF_rl_rt/",
  pattern = "\\.cf_rl_rt.nwk$|\\.nwk$|\\.treefile$|\\.tre$",
  full.names = TRUE
)

plot_list_node <- list()
stat_list_node <- list()  # store per-tree statistics rows

# length(tree_files)
# tree_files
# tr <- read.tree(tree_files[1])
# print(head(tr$node.label))

# parse taxonomic info from labels
parse_tip_info <- function(tip) {
  parts <- strsplit(tip, "~")[[1]]
  taxo  <- parts[2]
  taxo_split <- strsplit(taxo, "-")[[1]]
  
  superfamily <- taxo_split[1]
  family      <- ifelse(length(taxo_split) >= 2, taxo_split[2], NA)
  subfamily   <- ifelse(length(taxo_split) >= 3, taxo_split[3], NA)
  
  data.frame(
    tip = tip,
    superfamily = superfamily,
    family = family,
    subfamily = subfamily,
    stringsAsFactors = FALSE
  )
}

# extract gCF from label "UF/SH/gCF"
extract_gCF <- function(label) {
  if (is.na(label) || label == "") return(NA)
  
  parts <- strsplit(label, "/")[[1]]
  
  # concatenated tree: UFboot/SH-aLRT/gCF/sCF
  if (length(parts) == 4) {
    return(as.numeric(parts[3]))
  }
  
  # astral tree: LPP/gCF/sCF
  if (length(parts) == 3) {
    return(as.numeric(parts[2]))
  }
  
  return(NA)
}

# main loop
for (file in tree_files) {
  
  tr <- read.tree(file)
  if (is.null(tr$node.label)) next
  
  ntip <- length(tr$tip.label)
  tip_info <- do.call(rbind, lapply(tr$tip.label, parse_tip_info))
  
  # descendant tips
  desc_list <- vector("list", tr$Nnode)
  for (i in seq_len(tr$Nnode)) {
    node_id <- ntip + i
    desc_list[[i]] <- tr$tip.label[
      Descendants(tr, node_id, type = "tips")[[1]]
    ]
  }
  
  # build node df
  node_df <- data.frame(
    node = ntip + seq_len(tr$Nnode),
    gCF  = sapply(tr$node.label, extract_gCF),
    NodeType = NA_character_,
    IsStaphylinidae = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Classify internal nodes by taxonomic level.
  for (i in seq_len(tr$Nnode)) {
    
    tips <- desc_list[[i]]
    info <- tip_info[match(tips, tip_info$tip), ]
    
    sf <- unique(info$superfamily)
    fa <- unique(info$family)
    su <- unique(info$subfamily)
    
    n_sf <- length(sf[!is.na(sf)])
    n_fa <- length(fa[!is.na(fa)])
    n_su <- length(su[!is.na(su)])
    
    if (n_su == 1) {
      node_df$NodeType[i] <- "Within subfamily"
    } else if (n_fa == 1 && n_su > 1) {
      node_df$NodeType[i] <- "Subfamily"
    } else if (n_sf == 1 && n_fa > 1) {
      node_df$NodeType[i] <- "Family"
    }
    
    if (n_fa == 1 && fa[1] == "Staphylinidae") {
      node_df$IsStaphylinidae[i] <- TRUE
    }
  }
  
  node_df <- node_df %>%
    filter(!is.na(NodeType) & !is.na(gCF)) %>%
    mutate(
      NodeType = factor(
        NodeType,
        levels = c("Family", "Subfamily", "Within subfamily")
      )
    )
  
  if (nrow(node_df) == 0) next
  
  # Summarize gCF by hierarchical level.
  
  # Map plotting labels to table labels.
  node_df_stat <- node_df %>%
    mutate(
      LevelForTable = case_when(
        NodeType == "Family" ~ "family",
        NodeType == "Subfamily" ~ "subfamily",
        NodeType == "Within subfamily" ~ "within.sf",
        TRUE ~ NA_character_
      )
    )
  
  # Statistics for the three hierarchical levels.
  level_stat <- node_df_stat %>%
    filter(!is.na(LevelForTable)) %>%
    group_by(LevelForTable) %>%
    summarise(
      mean = mean(gCF),
      max  = max(gCF),
      min  = min(gCF),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = LevelForTable,
      values_from = c(mean, max, min),
      names_glue  = "gCF_{LevelForTable}_{.value}"
    )
  
  # Split subfamily-level nodes into Staphylinidae and non-Staphylinidae.
  sub_df <- node_df_stat %>% filter(NodeType == "Subfamily")
  
  sub_sta <- sub_df %>% filter(IsStaphylinidae)
  sub_non <- sub_df %>% filter(!IsStaphylinidae)
  
  sf_sta_mean <- if (nrow(sub_sta) > 0) mean(sub_sta$gCF) else NA_real_
  sf_sta_max  <- if (nrow(sub_sta) > 0) max(sub_sta$gCF)  else NA_real_
  sf_sta_min  <- if (nrow(sub_sta) > 0) min(sub_sta$gCF)  else NA_real_
  
  sf_non_mean <- if (nrow(sub_non) > 0) mean(sub_non$gCF) else NA_real_
  sf_non_max  <- if (nrow(sub_non) > 0) max(sub_non$gCF)  else NA_real_
  sf_non_min  <- if (nrow(sub_non) > 0) min(sub_non$gCF)  else NA_real_
  
  extra_stat <- data.frame(
    gCF_sf.Sta_mean = sf_sta_mean,
    gCF_sf.Sta_max  = sf_sta_max,
    gCF_sf.Sta_min  = sf_sta_min,
    gCF_sf.non.Sta_mean = sf_non_mean,
    gCF_sf.non.Sta_max  = sf_non_max,
    gCF_sf.non.Sta_min  = sf_non_min,
    stringsAsFactors = FALSE
  )
  
  stat_row <- bind_cols(
    data.frame(tree = basename(file), stringsAsFactors = FALSE),
    level_stat,
    extra_stat
  )
  
  stat_list_node[[basename(file)]] <- stat_row
  
  # Plot one tree.
  p <- ggplot(node_df, aes(x = NodeType, y = gCF)) +
    geom_jitter(
      aes(color = IsStaphylinidae),
      width = 0.2,
      height = 0,
      size = 2
    ) +
    scale_color_manual(
      values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
      guide = "none"
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      minor_breaks = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4)
    ) +
    labs(
      title = basename(file),
      x = "Hierarchical level",
      y = "gCF"
    )
  
  out_png <- paste0(
    "gCF.vs.node/",
    tools::file_path_sans_ext(basename(file)),
    "_gCF_vs_node.png"
  )
  
  ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  plot_list_node[[basename(file)]] <- p
}

# Write one summary row per tree.
if (length(stat_list_node) > 0) {
  
  gCF_stat_all <- bind_rows(stat_list_node)
  
  # Ensure all expected columns exist (fill with NA if missing)
  expected_cols <- c(
    "tree",
    "gCF_family_mean", "gCF_family_max", "gCF_family_min",
    "gCF_subfamily_mean", "gCF_subfamily_max", "gCF_subfamily_min",
    "gCF_within.sf_mean", "gCF_within.sf_max", "gCF_within.sf_min",
    "gCF_sf.Sta_mean", "gCF_sf.Sta_max", "gCF_sf.Sta_min",
    "gCF_sf.non.Sta_mean", "gCF_sf.non.Sta_max", "gCF_sf.non.Sta_min"
  )
  
  for (cc in expected_cols) {
    if (!cc %in% colnames(gCF_stat_all)) gCF_stat_all[[cc]] <- NA
  }
  
  gCF_stat_all <- gCF_stat_all %>%
    select(all_of(expected_cols))
  
  write.csv(
    gCF_stat_all,
    file = "gCF.vs.node/gCF_by.taxa.csv",
    row.names = FALSE
  )
}

# Combined panel.

dataset_names <- c(
  "ALL_matrix50","ALL_matrix50_len70","ALL_matrix75",
  "ALL_matrix75_len70","ALL_matrix90","ALL_matrix90_len70",
  "SRH_matrix50", "SRH_matrix50_len70","SRH_matrix75",
  "SRH_matrix75_len70","SRH_matrix90","SRH_matrix90_len70",
  "FLT_matrix50", "FLT_matrix75"
)

method_levels <- c("LG", "GHOST", "EXEHO", "PMSF", "Dayhoff6")

if (length(plot_list_node) > 0) {
  
  plot_info_node <- data.frame(
    file = names(plot_list_node),
    base = tools::file_path_sans_ext(names(plot_list_node)),
    stringsAsFactors = FALSE
  )
  
  plot_info_node$Method <- sub("_.*$", "", plot_info_node$base)
  plot_info_node$Dataset <- NA_character_
  
  for (ds in dataset_names) {
    hit <- grepl(ds, plot_info_node$base, fixed = TRUE)
    plot_info_node$Dataset[hit] <- ds
  }
  
  plot_info_node <- subset(
    plot_info_node,
    !is.na(Dataset) & Method %in% method_levels
  )
  
  plot_info_node$Method  <- factor(plot_info_node$Method,
                                   levels = method_levels)
  plot_info_node$Dataset <- factor(plot_info_node$Dataset,
                                   levels = dataset_names)
  
  row_panels <- list()
  
  for (ds in dataset_names) {
    sub_df <- subset(plot_info_node, Dataset == ds)
    
    row_plots <- vector("list", length(method_levels))
    names(row_plots) <- method_levels
    
    for (m in method_levels) {
      this <- sub_df[sub_df$Method == m, ]
      if (nrow(this) == 1) {
        row_plots[[m]] <- plot_list_node[[this$file]]
      } else {
        row_plots[[m]] <- ggplot() + theme_void()
      }
    }
    
    row_panels[[ds]] <- wrap_plots(
      row_plots,
      ncol = length(method_levels)
    )
  }
  
  final_panel_node <- wrap_plots(row_panels, ncol = 1)
  
  ggsave(
    filename = "gCF.vs.node/gCF_vs_node_panel.png",
    plot     = final_panel_node,
    width    = 30,
    height   = 4 * length(dataset_names),
    dpi      = 300,
    limitsize = FALSE
  )
}


# 5.5.3 sCF by hierarchical level
library(ape)
library(phangorn)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

dir.create("3_sCF.vs.node", showWarnings = FALSE)

tree_files <- list.files(
  "0_treefiles_CF_rl_rt/",
  pattern = "\\.nwk$|\\.treefile$|\\.tre$",
  full.names = TRUE
)

plot_list_node <- list()
stat_list_node <- list()   # NEW: store per-tree statistics

# parse taxonomic information from tip labels
parse_tip_info <- function(tip) {
  parts <- strsplit(tip, "~")[[1]]
  taxo  <- parts[2]
  taxo_split <- strsplit(taxo, "-")[[1]]
  
  superfamily <- taxo_split[1]
  family      <- ifelse(length(taxo_split) >= 2, taxo_split[2], NA)
  subfamily   <- ifelse(length(taxo_split) >= 3, taxo_split[3], NA)
  
  data.frame(
    tip = tip,
    superfamily = superfamily,
    family = family,
    subfamily = subfamily,
    stringsAsFactors = FALSE
  )
}

# extract sCF from label "UF/SH-aLRT/gCF/sCF"
extract_sCF <- function(label) {
  if (is.na(label) || label == "") return(NA)
  parts <- strsplit(label, "/")[[1]]
  if (length(parts) != 4) return(NA)
  as.numeric(parts[4])
}

# main loop
for (file in tree_files) {
  
  tr <- read.tree(file)
  if (is.null(tr$node.label)) next
  
  ntip <- length(tr$tip.label)
  tip_info <- do.call(rbind, lapply(tr$tip.label, parse_tip_info))
  
  # descendant tips
  desc_list <- vector("list", tr$Nnode)
  for (i in seq_len(tr$Nnode)) {
    node_id <- ntip + i
    desc_list[[i]] <- tr$tip.label[
      Descendants(tr, node_id, type = "tips")[[1]]
    ]
  }
  
  # build node-level data frame
  node_df <- data.frame(
    node = ntip + seq_len(tr$Nnode),
    sCF  = sapply(tr$node.label, extract_sCF),
    NodeType = NA_character_,
    IsStaphylinidae = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Classify internal nodes by taxonomic level.
  for (i in seq_len(tr$Nnode)) {
    
    tips <- desc_list[[i]]
    info <- tip_info[match(tips, tip_info$tip), ]
    
    sf <- unique(info$superfamily)
    fa <- unique(info$family)
    su <- unique(info$subfamily)
    
    n_sf <- length(sf[!is.na(sf)])
    n_fa <- length(fa[!is.na(fa)])
    n_su <- length(su[!is.na(su)])
    
    if (n_su == 1) {
      node_df$NodeType[i] <- "Within subfamily"
    } else if (n_fa == 1 && n_su > 1) {
      node_df$NodeType[i] <- "Subfamily"
    } else if (n_sf == 1 && n_fa > 1) {
      node_df$NodeType[i] <- "Family"
    }
    
    if (n_fa == 1 && fa[1] == "Staphylinidae") {
      node_df$IsStaphylinidae[i] <- TRUE
    }
  }
  
  node_df <- node_df %>%
    filter(!is.na(NodeType) & !is.na(sCF)) %>%
    mutate(
      NodeType = factor(
        NodeType,
        levels = c("Family", "Subfamily", "Within subfamily")
      ),
      Low_sCF = sCF < 33.33
    )
  
  if (nrow(node_df) == 0) next
  
  # Summarize sCF by hierarchical level.
  
  node_df_stat <- node_df %>%
    mutate(
      LevelForTable = case_when(
        NodeType == "Family" ~ "family",
        NodeType == "Subfamily" ~ "subfamily",
        NodeType == "Within subfamily" ~ "within.sf",
        TRUE ~ NA_character_
      )
    )
  
  # Statistics for the three hierarchical levels.
  level_stat <- node_df_stat %>%
    filter(!is.na(LevelForTable)) %>%
    group_by(LevelForTable) %>%
    summarise(
      mean = mean(sCF),
      max  = max(sCF),
      min  = min(sCF),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = LevelForTable,
      values_from = c(mean, max, min),
      names_glue  = "sCF_{LevelForTable}_{.value}"
    )
  
  # Split subfamily-level nodes into Staphylinidae and non-Staphylinidae.
  sub_df  <- node_df_stat %>% filter(NodeType == "Subfamily")
  sub_sta <- sub_df %>% filter(IsStaphylinidae)
  sub_non <- sub_df %>% filter(!IsStaphylinidae)
  
  extra_stat <- data.frame(
    sCF_sf.Sta_mean     = if (nrow(sub_sta) > 0) mean(sub_sta$sCF) else NA,
    sCF_sf.Sta_max      = if (nrow(sub_sta) > 0) max(sub_sta$sCF)  else NA,
    sCF_sf.Sta_min      = if (nrow(sub_sta) > 0) min(sub_sta$sCF)  else NA,
    sCF_sf.non.Sta_mean = if (nrow(sub_non) > 0) mean(sub_non$sCF) else NA,
    sCF_sf.non.Sta_max  = if (nrow(sub_non) > 0) max(sub_non$sCF)  else NA,
    sCF_sf.non.Sta_min  = if (nrow(sub_non) > 0) min(sub_non$sCF)  else NA,
    stringsAsFactors = FALSE
  )
  
  stat_row <- bind_cols(
    data.frame(tree = basename(file), stringsAsFactors = FALSE),
    level_stat,
    extra_stat
  )
  
  stat_list_node[[basename(file)]] <- stat_row
  
  # Plot one tree.
  p <- ggplot(node_df, aes(x = NodeType, y = sCF)) +
    geom_jitter(
      aes(color = IsStaphylinidae, alpha = Low_sCF),
      width = 0.2,
      height = 0,
      size = 2
    ) +
    geom_hline(
      yintercept = 33.33,
      linetype   = "dashed",
      color      = "darkgrey",
      linewidth  = 0.4
    ) +
    scale_color_manual(
      values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c("TRUE" = 0.5, "FALSE" = 1),
      guide = "none"
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 33.3, 50, 75, 100)
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = basename(file),
      x = "Hierarchical level",
      y = "sCF"
    )
  
  out_png <- paste0(
    "3_sCF.vs.node/",
    tools::file_path_sans_ext(basename(file)),
    "_sCF_vs_node.png"
  )
  
  ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  
  plot_list_node[[basename(file)]] <- p
}

# Write one summary row per tree.
if (length(stat_list_node) > 0) {
  
  sCF_stat_all <- bind_rows(stat_list_node)
  
  expected_cols <- c(
    "tree",
    "sCF_family_mean", "sCF_family_max", "sCF_family_min",
    "sCF_subfamily_mean", "sCF_subfamily_max", "sCF_subfamily_min",
    "sCF_within.sf_mean", "sCF_within.sf_max", "sCF_within.sf_min",
    "sCF_sf.Sta_mean", "sCF_sf.Sta_max", "sCF_sf.Sta_min",
    "sCF_sf.non.Sta_mean", "sCF_sf.non.Sta_max", "sCF_sf.non.Sta_min"
  )
  
  for (cc in expected_cols) {
    if (!cc %in% colnames(sCF_stat_all)) sCF_stat_all[[cc]] <- NA
  }
  
  sCF_stat_all <- sCF_stat_all %>%
    select(all_of(expected_cols))
  
  write.csv(
    sCF_stat_all,
    file = "3_sCF.vs.node/sCF_by.taxa.csv",
    row.names = FALSE
  )
}

# Combined panel.

dataset_names <- c(
  "ALL_matrix50","ALL_matrix50_len70","ALL_matrix75",
  "ALL_matrix75_len70","ALL_matrix90","ALL_matrix90_len70",
  "SRH_matrix50", "SRH_matrix50_len70","SRH_matrix75",
  "SRH_matrix75_len70","SRH_matrix90","SRH_matrix90_len70",
  "FLT_matrix50", "FLT_matrix75"
)

method_levels <- c("LG", "GHOST", "EXEHO", "PMSF", "Dayhoff6")

if (length(plot_list_node) > 0) {
  
  plot_info_node <- data.frame(
    file = names(plot_list_node),
    base = tools::file_path_sans_ext(names(plot_list_node)),
    stringsAsFactors = FALSE
  )
  
  plot_info_node$Method <- sub("_.*$", "", plot_info_node$base)
  plot_info_node$Dataset <- NA_character_
  
  for (ds in dataset_names) {
    hit <- grepl(ds, plot_info_node$base, fixed = TRUE)
    plot_info_node$Dataset[hit] <- ds
  }
  
  plot_info_node <- subset(
    plot_info_node,
    !is.na(Dataset) & Method %in% method_levels
  )
  
  plot_info_node$Method  <- factor(plot_info_node$Method,
                                   levels = method_levels)
  plot_info_node$Dataset <- factor(plot_info_node$Dataset,
                                   levels = dataset_names)
  
  row_panels <- list()
  
  for (ds in dataset_names) {
    
    sub_df <- subset(plot_info_node, Dataset == ds)
    
    row_plots <- vector("list", length(method_levels))
    names(row_plots) <- method_levels
    
    for (m in method_levels) {
      this <- sub_df[sub_df$Method == m, ]
      if (nrow(this) == 1) {
        row_plots[[m]] <- plot_list_node[[this$file]]
      } else {
        row_plots[[m]] <- ggplot() + theme_void()
      }
    }
    
    row_panels[[ds]] <- wrap_plots(
      row_plots,
      ncol = length(method_levels)
    )
  }
  
  final_panel_node <- wrap_plots(row_panels, ncol = 1)
  
  ggsave(
    filename = "3_sCF.vs.node/sCF_vs_node_panel.png",
    plot     = final_panel_node,
    width    = 30,
    height   = 4 * length(dataset_names),
    dpi      = 300,
    limitsize = FALSE
  )
}


# 5.5.4 Node-support trends across filtering strategies
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# 1. Create output folder
dir.create("5_node_support_concat", showWarnings = FALSE)

# 2. Read summary table
df <- read.csv("1_violin_CF/node_support_summary_CF.csv", stringsAsFactors = FALSE)

# 3. Keep PMSF trees only
df_pmsf <- df %>%
  filter(grepl("^PMSF_", Tree))

# 4. Function to parse x-axis group for a given taxa level
get_group <- function(tree_name, taxa_level) {
  
  if (grepl(paste0("_ALL_", taxa_level, "_len70_"), tree_name)) {
    return("ALL_len70")
  }
  if (grepl(paste0("_ALL_", taxa_level, "_"), tree_name) &&
      !grepl(paste0("_ALL_", taxa_level, "_len70_"), tree_name)) {
    return("ALL")
  }
  if (grepl(paste0("_SRH_", taxa_level, "_len70_"), tree_name)) {
    return("SRH_len70")
  }
  if (grepl(paste0("_SRH_", taxa_level, "_"), tree_name) &&
      !grepl(paste0("_SRH_", taxa_level, "_len70_"), tree_name)) {
    return("SRH")
  }
  if (grepl(paste0("_FLT_", taxa_level, "_"), tree_name)) {
    return("FLT")
  }
  
  return(NA)
}

# 5. Function to make one plot for one taxa level
make_node_support_plot <- function(data, taxa_level) {
  
  sub <- data %>%
    filter(grepl(taxa_level, Tree)) %>%
    mutate(group = sapply(Tree, get_group, taxa_level = taxa_level))
  
  x_order <- c("ALL", "ALL_len70", "SRH", "SRH_len70", "FLT")
  
  # Convert to long format
  plot_df <- bind_rows(
    sub %>%
      transmute(
        Tree,
        group,
        metric = "UFboot",
        mean = Mean_UFboot,
        min = Min_UFboot,
        max = Max_UFboot
      ),
    sub %>%
      transmute(
        Tree,
        group,
        metric = "SH_aLRT",
        mean = Mean_SHaLRT,
        min = Min_SHaLRT,
        max = Max_SHaLRT
      ),
    sub %>%
      transmute(
        Tree,
        group,
        metric = "gCF",
        mean = Mean_gCF,
        min = Min_gCF,
        max = Max_gCF
      ),
    sub %>%
      transmute(
        Tree,
        group,
        metric = "sCF",
        mean = Mean_sCF,
        min = Min_sCF,
        max = Max_sCF
      )
  ) %>%
    mutate(
      group = factor(group, levels = x_order),
      metric = factor(metric, levels = c("UFboot", "SH_aLRT", "gCF", "sCF"))
    )
  
  dodge_width <- 0.6
  
  p <- ggplot(
    plot_df,
    aes(
      x = group,
      y = mean,
      color = metric,
      shape = metric
    )
  ) +
    geom_errorbar(
      aes(ymin = min, ymax = max),
      position = position_dodge(width = dodge_width),
      width = 0.6,
      linewidth = 0.8
    ) +
    geom_point(
      position = position_dodge(width = dodge_width),
      size = 2.5
    ) +
    scale_color_manual(
      values = c(
        "UFboot"  = "firebrick",
        "SH_aLRT" = "gold3",
        "gCF"     = "#8499ce",
        "sCF"     = "#574774"
      ),
      breaks = c("UFboot", "SH_aLRT", "gCF", "sCF"),
      labels = c("UFboot", "SH-aLRT", "gCF", "sCF")
    ) +
    scale_shape_manual(
      values = c(
        "UFboot"  = 16,
        "SH_aLRT" = 17,
        "gCF"     = 18,
        "sCF"     = 15
      ),
      breaks = c("UFboot", "SH_aLRT", "gCF", "sCF"),
      labels = c("UFboot", "SH-aLRT", "gCF", "sCF")
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.title = element_blank(),
      aspect.ratio = 1
    ) +
    labs(
      title = paste("PMSF", taxa_level),
      x = "Filtering criteria / length coverage",
      y = "Node support"
    )
  
  return(p)
}

# 6. Make and save three plots
taxa_levels <- c("taxa50", "taxa75", "taxa90")

for (tx in taxa_levels) {
  
  p <- make_node_support_plot(df_pmsf, tx)
  
  ggsave(
    filename = paste0("5_node_support_concat/PMSF_", tx, "_node_support.png"),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
}

